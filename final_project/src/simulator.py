import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from src.state import NodeState
from src.reservoir import Reservoir
from src.well import Well
from src.pipe import Pipe
from src.compressor import DCS


class FieldSimulator:
    """
    Симулятор газового куста. Решает систему уравнений гидравлического равновесия
    и моделирует истощение пласта по времени.
    """

    def __init__(self, reservoir: Reservoir, wells: list, shlyf: Pipe, dcs: DCS):
        self.reservoir = reservoir
        self.wells = wells
        self.shlyf = shlyf
        self.dcs = dcs

    def _residuals(self, x: list, P_res: float) -> list:
        """
        Вектор невязок для fsolve.
        x = [q1, q2, q3, P_man]
        """
        q1, q2, q3, P_man = x
        qs = [q1, q2, q3]
        res = [0.0] * 4

        # === Уравнения 1-3: баланс для каждой скважины ===
        for i, well in enumerate(self.wells):
            # Итерационный поиск P_bhp для расчёта свойств в трубе
            P_bhp = P_man + 5.0  # начальное приближение: +5 атм к устьевому
            for _ in range(5):  # обычно хватает 2-3 итераций
                # Свойства при среднем давлении в трубе
                P_avg = (P_bhp + P_man) / 2.0
                dP_tube = well.pipe.dp(P_avg, qs[i]).dP
                P_bhp_new = P_man + dP_tube
                if abs(P_bhp_new - P_bhp) < 0.01:
                    break
                P_bhp = P_bhp_new

            # Коэффициент продуктивности при пластовом давлении
            C = well.get_productivity_index(P_res)

            # Невязка: q_i - C_i * (P_res - P_bhp) = 0
            # где P_bhp = P_man + dP_tube
            res[i] = qs[i] - C * (P_res - P_bhp)

        # === Уравнение 4: баланс шлейфа и ДКС ===
        q_total_sys = sum(qs) + self.dcs.q_ext
        # Для шлейфа: вход = P_man, выход = P_in_DCS
        dP_shlyf = self.shlyf.dp(P_man, q_total_sys).dP
        P_in_dcs = self.dcs.P_in()

        # Невязка: фактическое P_man должно равняться P_in_DCS + потери в шлейфе
        res[3] = P_man - (P_in_dcs + dP_shlyf)

        return res

    def solve(self, P_res: float) -> dict:
        """
        Находит рабочую точку системы при заданном пластовом давлении.
        Возвращает словарь NodeState для всех элементов.
        """
        # Начальное приближение по ТЗ
        x0 = [500.0, 500.0, 500.0, self.dcs.P_in() + 5.0]

        sol, info, ier, mesg = fsolve(self._residuals, x0, args=(P_res,), full_output=True)
        if ier != 1:
            print(f"[solve] Warning: fsolve convergence issue. {mesg}")

        q1, q2, q3, P_man = sol
        # Требование ТЗ: обнуляем отрицательные дебиты
        qs = [max(0.0, q) for q in [q1, q2, q3]]
        q_total_wells = sum(qs)

        states = {}
        for i, q in enumerate(qs):
            well = self.wells[i]
            # Пересчитываем точное состояние трубы с обнулённым дебитом
            P_in_tube_approx = P_man + 10.0
            pipe_state = well.pipe.dp(P_in_tube_approx, q)

            # Реальное забойное давление = устьевое + потери
            actual_P_bhp = P_man + pipe_state.dP

            states[f'well_{i + 1}'] = NodeState(
                name=f'well_{i + 1}',
                P_in=actual_P_bhp,
                P_out=P_man,
                dP=pipe_state.dP,
                q_std=q,
                q_res=pipe_state.q_res,
                v=pipe_state.v,
                rho=pipe_state.rho
            )

        # Состояние шлейфа
        q_total_sys = q_total_wells + self.dcs.q_ext
        states['shlyf'] = self.shlyf.dp(P_man, q_total_sys)

        # Состояние ДКС
        p_in_dcs = states['shlyf'].P_out
        p_out_dcs = self.dcs.P_line
        states['dcs'] = NodeState(
            name='dcs',
            P_in=p_in_dcs,
            P_out=p_out_dcs,
            dP=p_out_dcs - p_in_dcs,
            q_std=q_total_sys,
            q_res=None, v=None, rho=None
        )

        return states

    def run(self, N_days: int, dt: float = 1.0) -> pd.DataFrame:
        """
        Запускает динамическую симуляцию на N_days шагов.
        Возвращает DataFrame с историей разработки.
        """
        records = []
        current_Gp = 0.0

        for day in range(N_days):
            # Вывод прогресса каждые 10 шагов
            if day % 10 == 0 or day == N_days - 1:
                print(f"Day {day:3d}/{N_days} | P_res: {self.reservoir.resprops.P:.2f} атм")

            P_res = self.reservoir.resprops.P
            states = self.solve(P_res)

            q1 = states['well_1'].q_std
            q2 = states['well_2'].q_std
            q3 = states['well_3'].q_std
            P_man = states['well_1'].P_out  # У всех скважин общее устьевое/манифольдное давление

            q_total = q1 + q2 + q3
            current_Gp += q_total * dt

            records.append({
                't': day * dt,
                'P_res': P_res,
                'P_man': P_man,
                'q1': q1,
                'q2': q2,
                'q3': q3,
                'q_total': q_total,
                'Gp': current_Gp / 1000.0  # тыс. ст.м³
            })

            # Обновление пластового давления по материальному балансу
            # ВНИМАНИЕ: в баланс входит ТОЛЬКО дебит скважин. Сторонний газ q_ext не учитывается.
            P_next = self.reservoir.p2(q_total, dt)
            self.reservoir.resprops.P = P_next

        return pd.DataFrame(records)
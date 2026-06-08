import math
import numpy as np
from scipy.optimize import fsolve
from state import NodeState
from fluid import Fluid


class Pipe:
    """
    Гидравлическая модель трубы (НКТ, шлейф)
    Пошаговый расчёт с учётом изменения PVT-свойств вдоль ствола
    """

    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid,
                 vertical_depth: float = 0.0, name: str = ""):
        self.L = L
        self.D = D
        self.roughness = roughness
        self.fluid = fluid
        self.H = vertical_depth
        self.name = name
        self.g = 9.81

    def _calc_lambda(self, Re: float, rel_rough: float) -> float:
        """Расчёт коэффициента гидравлического сопротивления λ"""
        if Re < 1e-6:
            return 0.0
        if Re < 2300.0:
            return 64.0 / Re

        # Уравнение Колбрука-Уайта
        lam = 0.02
        for _ in range(50):
            term = rel_rough / 3.7 + 2.51 / (Re * math.sqrt(lam))
            lam_new = 1.0 / (-2.0 * math.log10(term)) ** 2
            if abs(lam_new - lam) < 1e-6:
                return lam_new
            lam = lam_new
        return lam

    def _calc_pressure_step(self, P: float, q: float, dl: float, dz: float) -> float:
        """
        Расчёт перепада давления на одном шаге (снизу вверх)
        Возвращает давление на следующем шаге
        """
        if P <= 0.1:
            return 0.1

        # PVT-свойства
        try:
            rho = self.fluid.ro(P)
            Bg = self.fluid.bg(P)
            mu = self.fluid.mu(P) / 1000.0  # сП -> Па·с

            # Защиты от некорректных значений
            rho = max(rho, 0.01)
            Bg = min(max(Bg, 1e-6), 1.0)
            mu = max(mu, 1e-6)
        except:
            return 0.1

        # Скорость потока
        A = math.pi * self.D ** 2 / 4.0
        v = (q / 86400.0) * Bg / A
        v = max(min(v, 300.0), 1e-5)

        # Число Рейнольдса
        Re = rho * v * self.D / mu
        Re = max(Re, 1.0)

        # Коэффициент трения
        rel_rough = self.roughness / self.D
        lam = self._calc_lambda(Re, rel_rough)
        lam = max(min(lam, 0.1), 0.008)

        # Перепады давления
        dp_fric = lam * (dl / self.D) * (rho * v ** 2 / 2.0)
        dp_grav = rho * self.g * dz
        dp_total = (dp_fric + dp_grav) / 101325.0  # Па -> атм

        # Защита от отрицательного давления
        P_new = P - dp_total
        return max(P_new, 0.1)

    def pwf_to_wh(self, P_bhp: float, q: float) -> float:
        """
        Расчёт устьевого давления по забойному (движение снизу вверх)
        """
        if P_bhp <= 0:
            return 0.1

        P = float(P_bhp)
        dl = 50.0  # УМЕНЬШИЛ шаг для плавности (было 100)
        n_steps = int(self.L // dl)
        if self.L % dl != 0:
            n_steps += 1

        for i in range(n_steps):
            # НЕ ПРЕРЫВАЕМ расчёт, даже если давление маленькое
            if P < 0.1:
                P = 0.1

            dl_i = dl if i != n_steps - 1 else self.L - dl * (n_steps - 1)
            dz = dl_i * (self.H / self.L) if self.L > 0 else 0

            P = self._calc_pressure_step(P, q, dl_i, dz)

        return max(P, 0.1)

    def wh_to_pwf(self, P_wh: float, q: float) -> float:
        """
        Расчёт забойного давления по устьевому (движение сверху вниз)
        Используется для метода dp
        """
        if P_wh <= 0:
            return 0.1

        P = float(P_wh)
        dl = 100.0
        n_steps = int(self.L // dl)
        if self.L % dl != 0:
            n_steps += 1

        for i in range(n_steps):
            if P <= 0.1:
                return 0.1

            dl_i = dl if i != n_steps - 1 else self.L - dl * (n_steps - 1)
            dz = dl_i * (self.H / self.L) if self.L > 0 else 0

            # Для движения сверху вниз добавляем перепад
            try:
                rho = self.fluid.ro(P)
                Bg = self.fluid.bg(P)
                mu = self.fluid.mu(P) / 1000.0

                rho = max(rho, 0.01)
                Bg = min(max(Bg, 1e-6), 1.0)
                mu = max(mu, 1e-6)
            except:
                return 0.1

            A = math.pi * self.D ** 2 / 4.0
            v = (q / 86400.0) * Bg / A
            v = max(min(v, 300.0), 1e-5)

            Re = rho * v * self.D / mu
            Re = max(Re, 1.0)

            rel_rough = self.roughness / self.D
            lam = self._calc_lambda(Re, rel_rough)
            lam = max(min(lam, 0.1), 0.008)

            dp_fric = lam * (dl_i / self.D) * (rho * v ** 2 / 2.0)
            dp_grav = rho * self.g * dz
            dp_total = (dp_fric + dp_grav) / 101325.0

            P = P + dp_total

        return max(P, 0.1)

    def dp(self, P_in: float, q: float) -> NodeState:
        """
        Основной метод. Расчёт перепада давления.
        Для НКТ: P_in = P_bhp (забой), P_out = P_wh (устье)
        Для шлейфа: P_in = P_out_shlyf (выход), P_out = P_man (вход)
        """
        if self.H > 0:
            # НКТ: движение снизу вверх (от забоя к устью)
            P_out = self.pwf_to_wh(P_in, q)
            dP = P_in - P_out
        else:
            # Шлейф: движение от выхода ко входу (против потока)
            # P_in здесь — это давление на выходе из шлейфа (P_out_shlyf)
            # Нам нужно найти давление на входе (P_man)
            P_out = self.wh_to_pwf(P_in, q)  # считаем от выхода ко входу
            dP = P_out - P_in  # P_man - P_out_shlyf

        # PVT-свойства при среднем давлении
        P_avg = (P_in + P_out) / 2.0
        try:
            rho = self.fluid.ro(P_avg)
            Bg = self.fluid.bg(P_avg)
            v = (4.0 * q * Bg) / (math.pi * self.D ** 2 * 86400.0)
        except:
            rho = 0.0
            Bg = 0.0
            v = 0.0

        return NodeState(
            name=self.name,
            P_in=P_in,
            P_out=P_out,
            dP=abs(dP),  # ВАЖНО: всегда положительный перепад
            q_std=q,
            q_res=q * Bg if Bg else 0.0,
            v=v,
            rho=rho
        )

    def get_vlp_point(self, P_wh: float, q: float) -> float:
        """
        Расчёт забойного давления при заданном устьевом и дебите.
        Используется для построения VLP.
        """
        return self.wh_to_pwf(P_wh, q)

    def get_vlp(self, P_man: float, q_values: list) -> tuple[list, list]:
        """
        Построение кривой VLP: P_bhp(q) при фиксированном P_man
        """
        qs = []
        pbhps = []
        for q in q_values:
            P_bhp = self.get_vlp_point(P_man, q)
            qs.append(q)
            pbhps.append(P_bhp)
        return qs, pbhps

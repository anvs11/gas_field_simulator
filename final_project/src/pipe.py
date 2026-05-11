import math
from src.state import NodeState
from src.fluid import Fluid


class Pipe:
    """
    Гидравлическая модель трубы (НКТ, шлейф)
    Расчёт перепада давления по уравнению Дарси–Вейсбаха
    """

    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid,
                 vertical_depth: float = 0.0, name: str = ""):
        """
        Параметры
        ----------
        L : float
            Длина трубы, м.
        D : float
            Внутренний диаметр, м.
        roughness : float
            Абсолютная шероховатость стенки δ, м.
        fluid : Fluid
            Объект с PVT-свойствами газа.
        vertical_depth : float
            Вертикальная глубина H, м (для горизонтального шлейфа = 0).
        name : str
            Идентификатор элемента (для NodeState).
        """
        self.L = L
        self.D = D
        self.roughness = roughness
        self.fluid = fluid
        self.H = vertical_depth
        self.name = name
        self.g = 9.81  # ускорение свободного падения, м/с²

    def _calc_lambda(self, Re: float, rel_rough: float) -> float:
        """
        Расчёт коэффициента гидравлического сопротивления λ
        """
        if Re < 2300.0: # при ламинарном режиме используется формула Пуазейля
            return 64.0 / Re

        # для турбулентного режима - неявное уравнение Колбрука–Уайта
        lam = 0.02  # начальное приближение
        for _ in range(50):  # с запасом, сойтись должно быстрее
            term = rel_rough / 3.7 + 2.51 / (Re * math.sqrt(lam))
            lam_new = 1.0 / (-2.0 * math.log10(term)) ** 2
            if abs(lam_new - lam) < 1e-6:
                return lam_new
            lam = lam_new

        return lam  # возврат последнего значения, если не сошлось за 50 шагов

    def dp(self, P: float, q: float) -> NodeState:
        """
        Основной метод. Вычисляет перепад давления и возвращает состояние элемента.

        Параметры
        ----------
        P : float
            Давление на входе в трубу, атм.
        q : float
            Коммерческий расход, ст.м³/сут.

        Возвращает
        -------
        NodeState
            Состояние потока с заполненными полями.
        """
        # свойства газа при давлении на входе
        rho = self.fluid.ro(P)  # плотность, кг/м³
        Bg = self.fluid.bg(P)  # объёмный коэффициент расширения газа, м³/ст.м³
        mu_cP = self.fluid.mu(P)  # сП
        mu_Pa_s = mu_cP / 1000.0  # Па·с (для числа Рейнольдса)

        # расчет скорости потока по формуле v = (4 * q_std * Bg) / (π * D² * 86400)
        v = (4.0 * q * Bg) / (math.pi * self.D ** 2 * 86400.0)  # м/с

        Re = (rho * v * self.D) / mu_Pa_s # число Рейнольдса

        # коэффициент трения
        rel_rough = self.roughness / self.D
        lam = self._calc_lambda(Re, rel_rough)

        # перепад давления (Па → атм)
        delta_P_friction_Pa = lam * (self.L / self.D) * (rho * v ** 2 / 2.0)
        delta_P_hydro_Pa = rho * self.g * self.H
        delta_P_total_Pa = delta_P_friction_Pa + delta_P_hydro_Pa
        delta_P_atm = delta_P_total_Pa / 101325.0

        P_in = P # давление на входе
        P_out = P_in - delta_P_atm # давление на выходе
        q_res = q * Bg  # объёмный расход при местных условиях

        return NodeState(
            name=self.name,
            P_in=P_in,
            P_out=P_out,
            dP=delta_P_atm,
            q_std=q,
            q_res=q_res,
            v=v,
            rho=rho
        )

    def get_vlp(self, P_man: float, q_values: list) -> tuple[list, list]:
        """
        Построение кривой VLP: P_bhp(q) при фиксированном давлении на манифольде

        Формула: P_bhp = P_man + ΔP_friction + ΔP_hydrostatic
        Свойства газа считаются при P_man (явная схема).

        Параметры
        ----------
        P_man : float
            Давление на манифольде/устье, атм.
        q_values : list
            Список дебитов для расчёта, ст.м³/сут.

        Возвращает
        ----------
        tuple[list, list]
            (q_std, P_bhp) — дебиты и соответствующие забойные давления
        """
        qs = []
        pbhps = []
        for q in q_values:
            # свойства при фиксированном P_man
            rho = self.fluid.ro(P_man)
            bg = self.fluid.bg(P_man)
            mu_Pa_s = self.fluid.mu(P_man) / 1000.0

            # расчет скорости
            v = (4.0 * q * bg) / (math.pi * self.D ** 2 * 86400.0)

            # определение Re и λ
            Re = (rho * v * self.D) / mu_Pa_s
            rel_rough = self.roughness / self.D
            lam = self._calc_lambda(Re, rel_rough)

            # перевод перепадов из Па в атм
            dP_fric = (lam * (self.L / self.D) * (rho * v ** 2 / 2.0)) / 101325.0
            dP_hydro = (rho * self.g * self.H) / 101325.0

            # забойное давление = устьевое + потери
            P_bhp = P_man + dP_fric + dP_hydro
            qs.append(q)
            pbhps.append(P_bhp)

        return qs, pbhps

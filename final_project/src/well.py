import math
from src.fluid import Fluid
from src.pipe import Pipe


class Well:
    """
    Модель скважины. Объединяет приток газа из пласта (закон Дарси)
    и подъём по НКТ (Pipe) для гидравлических расчётов
    """

    # Коэффициент перевода единиц в формуле Дарси для газа (мД, м, атм -> ст.м³/сут)
    BETA = 0.00852702

    def __init__(self, fluid: Fluid, k: float, h: float, re: float, rw: float, pipe: Pipe = None):
        """
        Параметры
        ----------
        fluid : Fluid
            Объект с PVT-свойствами.
        k : float
            Проницаемость, мД.
        h : float
            Эффективная мощность пласта, м.
        re : float
            Радиус контура питания, м.
        rw : float
            Радиус скважины, м.
        pipe : Pipe, optional
            Объект трубы (НКТ), по которой газ поднимается на поверхность.
        """
        self.fluid = fluid
        self.k = k
        self.h = h
        self.re = re
        self.rw = rw
        self.pipe = pipe

    def get_productivity_index(self, P_res: float) -> float:
        """
        Рассчитать коэффициент продуктивности C при текущем пластовом давлении.

        C = (BETA * k * h) / (mu(P_res) * ln(re / rw))

        Параметры
        ----------
        P_res : float
            Текущее пластовое давление, атм

        Возвращает
        ----------
        float
            Коэффициент продуктивности [ст.м³/(сут·атм)]
        """
        # Газ поступает из пласта, поэтому вязкость считается при пластовом давлении
        mu = self.fluid.mu(P_res)  # сП
        geo_term = self.k * self.h / math.log(self.re / self.rw) # геометрический фактор

        # C = BETA * (геометрия) / (вязкость)
        C = (self.BETA * geo_term) / mu # коэффициент продуктивности

        return C

    def q(self, P_res: float, P_bhp: float) -> float:
        """
        Рассчитать дебит скважины по закону Дарси (IPR).

        q_std = C * (P_res - P_bhp)

        Параметры
        ----------
        P_res : float
            Пластовое давление, атм.
        P_bhp : float
            Забойное давление, атм.

        Возвращает
        ----------
        float
            Дебит [ст.м³/сут].
        """
        # перестраховка от случая, если по ошибке забойное давление окажется выше пластового
        if P_bhp >= P_res:
            return 0.0

        C = self.get_productivity_index(P_res)
        flow_rate = C * (P_res - P_bhp)

        return flow_rate

    def get_ipr(self, P_res: float, pbhp_values: list) -> tuple[list, list]:
        """
        Вспомогательный метод для построения кривой притока (IPR)
        Рассчитывает дебиты для списка забойных давлений при фиксированном пластовом

        Параметры
        ----------
        P_res : float
            Фиксированное пластовое давление, атм
        pbhp_values : list
            Список значений P_bhp, атм

        Возвращает
        ----------
        tuple[list, list]
            (pbhp_values, q_values)
        """
        q_values = [self.q(P_res, pbhp) for pbhp in pbhp_values]
        return pbhp_values, q_values
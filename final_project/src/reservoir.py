from dataclasses import dataclass
from src.fluid import Fluid

@dataclass
class ResProps:
    """
    Контейнер параметров пласта
    Параметры продуктивности (k, h, re, rw) в Reservoir не хранятся
    """
    P: float # давление, атм
    V: float # объём, м³
    T: float # температура, К

class Reservoir:
    """
    "Бак", который хранит текущее состояние пласта ResProps
    и рассчитывает материальный баланс
    """
    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid = fluid

        # использую ту же логику расчёта, что и в Fluid
        P_std_pa = 101325.0  # стандартное давление, Па
        T_std = 293.15  # стандартная температура, K
        Z_std = 1.0  # коэффициент сверхжимаемости, безразмерный
        self.rho_std = (P_std_pa * fluid.M) / (Z_std * fluid.R * 1000.0 * T_std)

    def p2(self, q_total: float, dt: float = 1.0) -> float:
        """
        Рассчитать пластовое давление на следующем шаге по материальному балансу.

        Параметры
        ----------
        q_total : float
            Суммарный дебит всех скважин [ст.м³/сут].
        dt : float
            Шаг по времени [сут].

        Возвращает
        ----------
        float
            Новое давление пласта P_res [атм].

        Важно
        ----
        Метод не изменяет self.resprops.P. Обновление состояния
        выполняет FieldSimulator после получения результата.
        """
        # свойства газа при текущем пластовом давлении
        P_curr = self.resprops.P
        Z_curr = self.fluid.z(P_curr)
        rho_curr = self.fluid.ro(P_curr)

        # формула материального баланса (раздел 2.3)
        P_next = P_curr - (Z_curr * self.rho_std / rho_curr) * (q_total / self.resprops.V) * dt

        return P_next

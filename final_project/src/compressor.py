class DCS:
    def __init__(self, CR: float, P_line: float, q_ext: float = 0.0):
        """
        Параметры
        ----------
        CR : float
            Степень сжатия (>= 1.0).
            При CR = 1.0 возвращается P_line
        P_line : float
            Давление в магистральном газопроводе, атм.
        q_ext : float
            Расход стороннего газа, поступающего на манифолд, ст.м³/сут.
        """
        if CR < 1.0:
            raise ValueError("Степень сжатия CR должна быть >= 1.0")

        self.CR = CR
        self.P_line = P_line
        self.q_ext = q_ext

    def P_in(self) -> float:
        """
        Рассчитать давление на входе в ДКС.

        Формула: P_in = P_line / CR
        При CR = 1.0 ДКС отключена -> P_in = P_line

        Возвращает
        ----------
        float
            Давление на входе, атм.
        """
        if self.CR == 1.0:
            return self.P_line
        return self.P_line / self.CR
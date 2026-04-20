class LinearInterpolator:
    """
    Линейный интерполятор на чистом Python.

    Параметры
    ----------
    xs : list
        Узловые точки по X (отсортированы по возрастанию).
    ys : list
        Значения функции в узловых точках.
    """

    def __init__(self, xs: list, ys: list):
        if len(xs) != len(ys):
            raise ValueError(f"Длины xs и ys должны совпадать: {len(xs)} != {len(ys)}")

        if len(xs) < 2:
            raise ValueError("Для интерполяции нужно минимум 2 точки")

        # Проверка сортировки по возрастанию
        for i in range(len(xs) - 1):
            if xs[i] >= xs[i + 1]:
                raise ValueError(f"xs должны быть строго отсортированы по возрастанию")

        self.xs = xs
        self.ys = ys

    def predict(self, xp: float) -> float:
        """
        Вычислить интерполированное значение для заданного xp.

        Параметры
        ----------
        xp : float
            Точка, в которой нужно найти значение.

        Возвращает
        ----------
        float
            Интерполированное значение.

        Raises
        ------
        ValueError
            Если xp выходит за границы [xs[0], xs[-1]].
        """
        # Проверка границ
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError(
                f"Точка xp={xp} выходит за границы интерполяции "
                f"[{self.xs[0]}, {self.xs[-1]}]"
            )

        # Случай, когда xp точно совпадает с последней точкой
        if xp == self.xs[-1]:
            return self.ys[-1]

        # Поиск интервала, содержащего xp
        for i in range(len(self.xs) - 1):
            if self.xs[i] <= xp < self.xs[i + 1]:
                # Линейная интерполяция
                x0, x1 = self.xs[i], self.xs[i + 1]
                y0, y1 = self.ys[i], self.ys[i + 1]
                yp = y0 + (y1 - y0) * (xp - x0) / (x1 - x0)
                return yp
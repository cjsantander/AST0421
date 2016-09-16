import numpy as np


class TransParameters:

    def __init__(self, period, transit_depth, time_t, time_f):
        # Creamos los atributos.
        # Parámetros observables.
        self.p = period
        self.t_d = transit_depth
        self.t_t = time_t
        self.t_f = time_f

        # Cálculos comunes.
        self.sin_tf = np.power(np.sin(self.t_f * np.pi / self.p), 2)
        self.sin_tt = np.power(np.sin(self.t_t * np.pi / self.p), 2)
        self.sq_td = np.sqrt(self.t_d)

        # Parámetros saqueados.
        self.radius_planet = self.sq_td
        self.density_star = self._eq9()
        self.inclination = self._eq13()
        self.a_r = self._eq8()
        self.b = self._eq7()

    # Para obtener b.
    # def _eq7(self):
    #     sin = self.sin_tf / self.sin_tt

    #     num = np.power(1 - self.sq_td, 2)
    #     num -= np.power(1 + self.sq_td, 2) * sin

    #     den = 1 - sin

    #     return np.sqrt(num / den)
    def _eq7(self):
        cuo = np.power(self.t_f / self.t_t, 2)

        num = np.power(1 - self.sq_td, 2)
        num -= cuo * np.power(1 + self.sq_td, 2)

        den = 1 - cuo

        return np.sqrt(num / den)

    # Para obtener a/R(estrella)
    # def _eq8(self):
    #     b2 = np.power(self._eq7(), 2)

    #     num = np.power(1 + self.sq_td, 2)
    #     num -= b2 * (1 - self.sin_tt)

    #     den = self.sin_tt

    #     return np.sqrt(num / den)
    def _eq8(self):
        num = 2 * self.p * np.sqrt(self.sq_td)
        den = np.pi * np.sqrt(np.power(self.t_t, 2) - np.power(self.t_f, 2))

        return num / den

    # Para obtener la densidad(estrella)
    # def _eq9(self):
    #     const = np.power(365.25, 2) / np.power(215, 3)  # En base al Sol
    #     return const * np.power(self._eq8(), 3)
    def _eq9(self):
        num = np.sqrt(self.sq_td)
        num = np.power(num, 3)
        num = self.p * num

        den = np.power(self.t_t, 2) - np.power(self.t_f, 2)
        den = np.power(den, 3)
        den = np.sqrt(den)

        return 3.46 * np.power(10., -3) * num / den

    # Para obtener i
    def _eq13(self):
        arg = self._eq7() * np.power(self._eq8(), -1)
        return np.arccos(arg)

from __future__ import division
from .fourier import fourier_transform, inverse_fourier_transform
from .sfwexpr import sfwExpr
from .sym import fsym, ssym, tsym, omegasym, j


class omegaExpr(sfwExpr):

    """Fourier domain expression or symbol (angular frequency)."""

    var = omegasym
    domain_name = 'Angular frequency'
    domain_units = 'rad/s'

    def __init__(self, val, **assumptions):

        assumptions['real'] = True
        super(omegaExpr, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = tExpr

        if self.expr.find(ssym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find(tsym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on t' % self.expr)

    def inverse_fourier(self):
        """Attempt inverse Fourier transform."""

        from .symbols import f, pi        

        return self(2 * pi * f).inverse_fourier()

    def time(self):
        """Alias for inverse_fourier."""

        return self.inverse_fourier()

    def plot(self, wvector=None, **kwargs):
        """Plot angular frequency response at values specified by wvector.

        There are many plotting options, see matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_frequency=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.
        """

        from .plot import plot_angular_frequency
        return plot_angular_frequency(self, wvector, **kwargs)


class Yomega(omegaExpr):

    """omega-domain admittance."""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, **assumptions):

        super(Yomega, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Yt

    def cpt(self):
        from .oneport import G, C, L, Y

        if self.is_number:
            return G(self.expr)

        y = self * j * omega

        if y.is_number:
            return L((1 / y).expr)

        y = self / (j * omega)

        if y.is_number:
            return C(y.expr)

        return Y(self)


class Zomega(omegaExpr):

    """omega-domain impedance."""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, **assumptions):

        super(Zomega, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Zt

    def cpt(self):
        from .oneport import R, C, L, Z

        if self.is_number:
            return R(self.expr)

        z = self * j * omega

        if z.is_number:
            return C((1 / z).expr)

        z = self / (j * omega)

        if z.is_number:
            return L(z.expr)

        return Z(self)


class Vomega(omegaExpr):

    """omega-domain voltage (units V/rad/s)."""

    quantity = 'Voltage spectrum'
    units = 'V/rad/s'

    def __init__(self, val, **assumptions):

        super(Vomega, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Vt


class Iomega(omegaExpr):

    """omega-domain current (units A/rad/s)."""

    quantity = 'Current spectrum'
    units = 'A/rad/s'

    def __init__(self, val, **assumptions):

        super(Iomega, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = It


class Homega(omegaExpr):

    """omega-domain transfer function response."""

    quantity = 'Transfer function'
    units = ''

    def __init__(self, val, **assumptions):

        super(Homega, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Ht

from .texpr import Ht, It, Vt, Yt, Zt, tExpr
omega = omegaExpr('omega')


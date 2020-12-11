"""This module provides the omegaExpr class to represent omega-domain
(angular frequency Fourier domain) expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .fourier import inverse_fourier_transform
from .expr import Expr, expr
from .sym import fsym, ssym, tsym, omegasym, omega0sym, j, pi
from sympy import Expr as symExpr


class AngularFourierDomainExpression(Expr):

    """Fourier domain expression or symbol (angular frequency)."""

    var = omegasym
    domain_name = 'Angular frequency'
    domain_units = 'rad/s'

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainExpression, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainExpression

        if self.expr.find(ssym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find(tsym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on t' % self.expr)

    def inverse_fourier(self):
        """Attempt inverse Fourier transform."""

        expr = self.subs(2 * pi * fsym)
        result = inverse_fourier_transform(expr, fsym, tsym)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        else:
            result = TimeDomainExpression(result)
        return result

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

    def transform(self, arg, **assumptions):
        """Transform into a different domain."""

        from .fexpr import FourierDomainExpression, f

        arg = expr(arg)        
        if isinstance(arg, FourierDomainExpression):
            result = self.subs(2 * pi * f)
            return result.subs(arg, **assumptions)
        return super(AngularFourierDomainExpression, self).transform(arg, **assumptions)
        

class AngularFourierDomainAdmittance(AngularFourierDomainExpression):

    """omega-domain admittance."""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainAdmittance, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainAdmittance

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


class AngularFourierDomainImpedance(AngularFourierDomainExpression):

    """omega-domain impedance."""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainImpedance, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainImpedance

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


class AngularFourierDomainVoltage(AngularFourierDomainExpression):

    """omega-domain voltage (units V/rad/s)."""

    quantity = 'Voltage spectrum'
    units = 'V/rad/s'
    superkind = 'Voltage'    

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainVoltage, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainVoltage

    def __mul__(self, x):
        """Multiply"""

        if isinstance(x, AngularFourierDomainAdmittance):
            return AngularFourierDomainCurrent(super(AngularFourierDomainVoltage, self).__mul__(x))
        if isinstance(x, (ConstantExpression, AngularFourierDomainExpression, symExpr, int, float, complex)):
            return super(AngularFourierDomainVoltage, self).__mul__(x)
        self._incompatible(x, '*')

    def __truediv__(self, x):
        """Divide"""

        if isinstance(x, AngularFourierDomainImpedance):
            return AngularFourierDomainCurrent(super(AngularFourierDomainVoltage, self).__truediv__(x))
        if isinstance(x, (ConstantExpression, AngularFourierDomainExpression, symExpr, int, float, complex)):
            return super(AngularFourierDomainVoltage, self).__truediv__(x)
        self._incompatible(x, '/')
        

class AngularFourierDomainCurrent(AngularFourierDomainExpression):

    """omega-domain current (units A/rad/s)."""

    quantity = 'Current spectrum'
    units = 'A/rad/s'
    superkind = 'Current'        

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainCurrent, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainCurrent

    def __mul__(self, x):
        """Multiply"""

        if isinstance(x, AngularFourierDomainImpedance):
            return AngularFourierDomainVoltage(super(AngularFourierDomainCurrent, self).__mul__(x))            
        if isinstance(x, (ConstantExpression, AngularFourierDomainExpression, symExpr, int, float, complex)):
            return super(AngularFourierDomainCurrent, self).__mul__(x)
        self._incompatible(x, '*')        

    def __truediv__(self, x):
        """Divide"""

        if isinstance(x, AngularFourierDomainAdmittance):
            return AngularFourierDomainVoltage(super(AngularFourierDomainCurrent, self).__truediv__(x))
        if isinstance(x, (ConstantExpression, AngularFourierDomainExpression, symExpr, int, float, complex)):
            return super(AngularFourierDomainCurrent, self).__truediv__(x)
        self._incompatible(x, '/')                        

        
class AngularFourierDomainTransferFunction(AngularFourierDomainExpression):

    """omega-domain transfer function response."""

    quantity = 'Transfer function'
    units = ''

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainTransferFunction, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainImpulseResponse

        
def omegaexpr(arg, **assumptions):
    """Create omegaExpr object.  If `arg` is omegasym return omega"""

    if arg is omegasym:
        return omega
    if arg is omega0sym:
        return omega0    
    return AngularFourierDomainExpression(arg, **assumptions)

        
from .texpr import TimeDomainImpulseResponse, TimeDomainCurrent, TimeDomainVoltage, TimeDomainAdmittance, TimeDomainImpedance, TimeDomainExpression
from .cexpr import ConstantExpression
omega = AngularFourierDomainExpression('omega')
omega0 = AngularFourierDomainExpression('omega_0')

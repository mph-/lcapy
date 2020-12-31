"""This module provides the AngularFourierDomainExpression class to
represent omega-domain (angular frequency Fourier domain) expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .fourier import inverse_fourier_transform
from .expr import Expr, expr
from .sym import fsym, ssym, tsym, omegasym, omega0sym, j, pi
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin
from sympy import Expr as symExpr

__all__ = ('omegaexpr', )

class AngularFourierDomainExpression(Expr):
    """Fourier domain expression or symbol (angular frequency)."""

    var = omegasym
    domain = 'Angular Fourier'
    domain_label = 'Angular frequency'
    domain_units = 'rad/s'
    is_angular_fourier_domain = True
    is_transform_domain = True    

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainExpression, self).__init__(val, **assumptions)

        if self.expr.find(ssym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find(tsym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on t' % self.expr)

    def _class_by_quantity(self, quantity):

        if quantity == 'voltage':
            return AngularFourierDomainVoltage
        elif quantity == 'current':
            return AngularFourierDomainCurrent
        elif quantity == 'impedance':
            return AngularFourierDomainImpedance
        elif quantity == 'admittance':
            return AngularFourierDomainAdmittance
        elif quantity == 'transfer':
            return AngularFourierDomainTransferFunction
        elif quantity == 'undefined':
            return AngularFourierDomainExpression                                
        raise ValueError('Unknown quantity %s' % quantity)

    def as_expr(self):
        return AngularFourierDomainExpression(self)

    def as_voltage(self):
        return AngularFourierDomainVoltage(self)

    def as_current(self):
        return AngularFourierDomainCurrent(self)    

    def as_impedance(self):
        return AngularFourierDomainImpedance(self)

    def as_admittance(self):
        return AngularFourierDomainAdmittance(self)

    def as_transfer(self):
        return AngularFourierDomainTransferFunction(self)    
    
    def inverse_fourier(self, **assumptions):
        """Attempt inverse Fourier transform."""

        expr = self.subs(2 * pi * fsym)
        result = inverse_fourier_transform(expr, fsym, tsym)

        return self.wrap(TimeDomainExpression(result, **assumptions))                

    def time(self, **assumptions):
        """Alias for inverse_fourier."""

        return self.inverse_fourier(**assumptions)

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""

        return self
    
    def fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        
        from .symbols import f
        
        result = self.subs(2 * pi * f)
        return self.wrap(result)                

    def laplace(self, **assumptions):
        """Convert to Laplace domain."""

        result = self.time().laplace()
        return self.wrap(LaplaceDomainExpression(result, **assumptions))

    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        return self.time(**assumptions).phasor(**assumptions)

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


class AngularFourierDomainAdmittance(AdmittanceMixin, AngularFourierDomainExpression):
    """omega-domain admittance."""
    pass


class AngularFourierDomainImpedance(ImpedanceMixin, AngularFourierDomainExpression):
    """omega-domain impedance."""
    pass


class AngularFourierDomainVoltage(VoltageMixin, AngularFourierDomainExpression):
    """omega-domain voltage (units V/rad/s)."""

    quantity_label = 'Voltage spectrum'
    units = 'V/rad/s'


class AngularFourierDomainCurrent(CurrentMixin, AngularFourierDomainExpression):
    """omega-domain current (units A/rad/s)."""

    quantity_label = 'Current spectrum'
    units = 'A/rad/s'

        
class AngularFourierDomainTransferFunction(TransferMixin, AngularFourierDomainExpression):
    """omega-domain transfer function response."""
    pass

        
def omegaexpr(arg, **assumptions):
    """Create AngularFourierDomainExpression object.  If `arg` is omegasym return omega"""

    if arg is omegasym:
        return omega
    if arg is omega0sym:
        return omega0    
    return AngularFourierDomainExpression(arg, **assumptions)

        
from .texpr import TimeDomainExpression
from .sexpr import LaplaceDomainExpression
from .cexpr import ConstantExpression
omega = AngularFourierDomainExpression('omega')
omega0 = AngularFourierDomainExpression('omega_0')

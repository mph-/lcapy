"""This module provides the ConstantExpression class to represent constant expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from .expr import Expr
from .sym import symbols_find
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin

class ConstantExpression(Expr):
    """Constant real expression or symbol.

    If symbols in the expression are known to be negative, use
    ConstantExpression(expr, positive=False)

    """

    domain = 'Constant'
    is_constant_domain = True

    def __init__(self, val, **assumptions):

        symbols = symbols_find(val)
        for symbol in ('s', 'omega', 't', 'f'):
            if symbol in symbols:
                raise ValueError(
                    'constant expression %s cannot depend on %s' % (val, symbol))

        assumptions['dc'] = True        
        super(ConstantExpression, self).__init__(val, **assumptions)

    def _class_by_quantity(self, quantity):

        if quantity == 'voltage':
            return ConstantVoltage
        elif quantity == 'current':
            return ConstantCurrent
        elif quantity == 'impedance':
            return ConstantImpedance
        elif quantity == 'admittance':
            return ConstantAdmittance
        elif quantity == 'transfer':
            return ConstantTransferFunction
        elif quantity == 'undefined':
            return ConstantExpression                        
        raise ValueError('Unknown quantity %s' % quantity)
        
    def as_expr(self):
        return ConstantExpression(self)
        
    def as_voltage(self):
        return ConstantVoltage(self)

    def as_current(self):
        return ConstantCurrent(self)    

    def as_impedance(self):
        return ConstantImpedance(self)

    def as_admittance(self):
        return ConstantAdmittance(self)

    def as_transfer(self):
        return ConstantTransferFunction(self)
    
    def rms(self):
        return {ConstantVoltage: TimeDomainVoltage, ConstantCurrent : TimeDomainCurrent}[self.__class__](self)

    def fourier(self):
        """Convert to Fourier domain representation."""

        return self.time().fourier()

    def angular_fourier(self):
        """Convert to angular Fourier domain representation."""

        return self.time().angular_fourier()        

    def canonical(self, factor_const=True):
        # Minor optimisation
        return self

    def phasor(self, **assumptions):

        return self.laplace(**assumptions).phasor(**assumptions)
    
    
class ConstantVoltage(VoltageMixin, ConstantExpression):
    """This is considered a constant time-domain expression."""

    def cpt(self):
        from .oneport import Vdc
        return Vdc(self)

    def time(self, **assumptions):
        return TimeDomainVoltage(self, **assumptions)

    def laplace(self):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    
class ConstantCurrent(CurrentMixin, ConstantExpression):
    """This is considered a constant time-domain expression."""    

    def cpt(self):
        from .oneport import Idc
        return Idc(self)

    def time(self, **assumptions):
        return TimeDomainCurrent(self, **assumptions)

    def laplace(self):
        """Convert to Laplace domain representation."""

        return self.time().laplace()    

    
class ConstantImpedance(ImpedanceMixin, ConstantExpression):
    """This is considered a constant Laplace-domain expression."""

    def laplace(self, **assumptions):
        return LaplaceDomainImpedance(self, **assumptions)    

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.laplace(**assumptions).time(**assumptions)

    
class ConstantAdmittance(AdmittanceMixin, ConstantExpression):
    """This is considered a constant Laplace-domain expression."""    

    def laplace(self, **assumptions):
        return LaplaceDomainAdmittance(self, **assumptions)

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.laplace(**assumptions).time(**assumptions)    


class ConstantTransferFunction(TransferMixin, ConstantExpression):
    """This is considered a constant Laplace-domain expression."""    

    def laplace(self, **assumptions):
        return LaplaceDomainTransferFunction(self, **assumptions)

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.laplace(**assumptions).time(**assumptions)    


def cexpr(arg, **assumptions):
    """Create Lcapy constant expression from arg."""

    return ConstantExpression(arg)


from .texpr import t, TimeDomainCurrent, TimeDomainVoltage, TimeDomainExpression
from .sexpr import s, LaplaceDomainImpedance, LaplaceDomainAdmittance
from .sexpr import LaplaceDomainTransferFunction

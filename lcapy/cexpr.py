"""This module provides the ConstantExpression class to represent constant expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from .expr import Expr
from .sym import symbols_find
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin


class ConstantExpression(Expr):
    """Constant real expression or symbol.

    If symbols in the expression are known to be negative, use
    ConstantExpression(expr, positive=False)

    """

    domain = 'Constant'
    is_constant = True

    def __init__(self, val, **assumptions):

        symbols = symbols_find(val)
        for symbol in ('s', 'omega', 't', 'f'):
            if symbol in symbols:
                raise ValueError(
                    'constant expression %s cannot depend on %s' % (val, symbol))

        assumptions['dc'] = True        
        super(ConstantExpression, self).__init__(val, **assumptions)

    @classmethod
    def as_voltage(cls, expr, **assumptions):
        return ConstantVoltage(expr, **assumptions)

    @classmethod
    def as_current(cls, expr, **assumptions):
        return ConstantCurrent(expr, **assumptions)    

    # @classmethod
    # def as_impedance(cls, expr, **assumptions):
    #     return ConstantImpedance(expr)

    # @classmethod
    # def as_admittance(cls, expr):
    #     return ConstantAdmittance(expr)

    # @classmethod
    # def as_transfer(cls, expr):
    #     return ConstantTransferFunction(expr)

    def rms(self):
        return {ConstantVoltage: TimeDomainVoltage, ConstantCurrent : TimeDomainCurrent}[self.__class__](self)

    def laplace(self):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    def fourier(self):
        """Convert to Fourier domain representation."""

        return self.time().fourier()    

    def canonical(self, factor_const=True):
        # Minor optimisation
        return self

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.wrap(TimeDomainExpression(result, **assumptions))

    
class ConstantVoltage(VoltageMixin, ConstantExpression):

    def cpt(self):
        from .oneport import Vdc
        return Vdc(self)

    def time(self, **assumptions):
        return TimeDomainVoltage(self)

    
class ConstantCurrent(CurrentMixin, ConstantExpression):

    def cpt(self):
        from .oneport import Idc
        return Idc(self)

    def time(self, **assumptions):
        return TimeDomainCurrent(self)


from .texpr import t, TimeDomainCurrent, TimeDomainVoltage, TimeDomainExpression
from .sexpr import s

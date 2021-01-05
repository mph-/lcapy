"""This module provides the ConstantExpression class to represent constant expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from .expr import Expr
from .domains import ConstantDomain
from .sym import symbols_find
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin

class ConstantExpression(ConstantDomain, Expr):
    """Constant real expression or symbol.

    If symbols in the expression are known to be negative, use
    ConstantExpression(expr, positive=False)

    """

    def __init__(self, val, **assumptions):

        symbols = symbols_find(val)
        for symbol in ('s', 'omega', 't', 'f'):
            if symbol in symbols:
                raise ValueError(
                    'constant expression %s cannot depend on %s' % (val, symbol))

        assumptions['dc'] = True        
        super(ConstantExpression, self).__init__(val, **assumptions)

    def as_expr(self):
        return ConstantExpression(self)
        
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

    def time(self, **assumptions):
        return self.change(self, domain='time', **assumptions)

    def laplace(self):
        """Convert to Laplace domain representation."""

        return self.time().laplace()    
    

class ConstantTimeExpression(ConstantExpression):
    pass
    

class ConstantFrequencyExpression(ConstantExpression):

    def laplace(self, **assumptions):
        return self.change(self, domain='laplace', **assumptions)

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.laplace(**assumptions).time(**assumptions)        

    
class ConstantVoltage(VoltageMixin, ConstantTimeExpression):
    """This is considered a constant time-domain expression."""

    def cpt(self):
        from .oneport import Vdc
        return Vdc(self)

    
class ConstantCurrent(CurrentMixin, ConstantTimeExpression):
    """This is considered a constant time-domain expression."""    

    def cpt(self):
        from .oneport import Idc
        return Idc(self)

    
class ConstantImpedance(ImpedanceMixin, ConstantFrequencyExpression):
    """This is considered a constant Laplace-domain expression."""
    pass

    
class ConstantAdmittance(AdmittanceMixin, ConstantFrequencyExpression):
    """This is considered a constant Laplace-domain expression."""    
    pass


class ConstantTransferFunction(TransferMixin, ConstantFrequencyExpression):
    """This is considered a constant Laplace-domain expression."""    
    pass


def cexpr(arg, **assumptions):
    """Create Lcapy constant expression from arg."""

    return ConstantExpression(arg)


from .expressionclasses import expressionclasses

classes = expressionclasses.make(ConstantExpression)

classes['voltage'] = ConstantVoltage
classes['current'] = ConstantCurrent
classes['admittance'] = ConstantAdmittance
classes['impedance'] = ConstantImpedance
classes['transfer'] = ConstantTransferFunction
expressionclasses.add('constant', classes)

from .texpr import t, TimeDomainCurrent, TimeDomainVoltage, TimeDomainExpression
from .sexpr import s, LaplaceDomainExpression

"""This module provides the ConstantDomainExpression class to represent constant expressions.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from .expr import Expr, expr_make
from .domains import ConstantDomain
from .sym import symbols_find
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin

class ConstantDomainExpression(ConstantDomain, Expr):
    """Constant real expression or symbol.

    If symbols in the expression are known to be negative, use
    ConstantDomainExpression(expr, positive=False)

    """

    def __init__(self, val, **assumptions):

        symbols = symbols_find(val)
        for symbol in ('s', 'omega', 't', 'f'):
            if symbol in symbols:
                raise ValueError(
                    'constant expression %s cannot depend on %s' % (val, symbol))

        assumptions['dc'] = True        
        super(ConstantDomainExpression, self).__init__(val, **assumptions)

    def as_expr(self):
        return ConstantDomainExpression(self)

    def rms(self):
        """Return root mean square."""
        return self

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
    

class ConstantTimeDomainExpression(ConstantDomainExpression):

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined':
            return ConstantTimeDomainExpression

        return super(ConstantTimeDomainExpression, self)._class_by_quantity(quantity, domain)
    

class ConstantFrequencyDomainExpression(ConstantDomainExpression):

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined':
            return ConstantFrequencyDomainExpression

        return super(ConstantFrequencyDomainExpression, self)._class_by_quantity(quantity, domain)
    
    def laplace(self, **assumptions):
        return self.change(self, domain='laplace', **assumptions)

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.laplace(**assumptions).time(**assumptions)        

    
def cexpr(arg, **assumptions):
    """Create Lcapy constant expression from `arg`.

    By default, `arg` is assumed to be positive.  If symbols in the
    `arg` are known to be negative, use `cexpr(arg, positive=False)`.

    """
    try:
        quantity = arg.quantity
    except:
        quantity = 'undefined'

    if quantity == 'undefined':
        # Sit on the fence rather than choosing ConstantTimeDomainExpression
        return ConstantDomainExpression(arg, **assumptions)
    
    return expr_make('constant', arg, **assumptions)    


from .expressionclasses import expressionclasses

expressionclasses.register('constant', ConstantTimeDomainExpression, ConstantFrequencyDomainExpression,
                           ('voltage', 'current', 'voltagesquared', 'currentsquared', 'undefined'))

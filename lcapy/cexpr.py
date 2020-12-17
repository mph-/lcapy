from .expr import Expr
from .sym import symbols_find


class ConstantExpression(Expr):
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

    def rms(self):
        return {ConstantVoltage: TimeDomainVoltage, ConstantCurrent : TimeDomainCurrent}[self.__class__](self)

    def laplace(self):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    def canonical(self, factor_const=True):
        # Minor optimisation
        return self

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return TimeDomainExpression(self)

    
class ConstantVoltage(ConstantExpression):

    def __init__(self, val, **assumptions):

        super(ConstantVoltage, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = TimeDomainVoltage

    def cpt(self):
        from .oneport import Vdc
        return Vdc(self)

    def time(self, **assumptions):
        return TimeDomainVoltage(self)

    
class ConstantCurrent(ConstantExpression):

    def __init__(self, val, **assumptions):

        super(ConstantCurrent, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = TimeDomainCurrent

    def cpt(self):
        from .oneport import Idc
        return Idc(self)

    def time(self, **assumptions):
        return TimeDomainCurrent(self)


from .texpr import t, TimeDomainCurrent, TimeDomainVoltage, TimeDomainExpression
from .sexpr import s

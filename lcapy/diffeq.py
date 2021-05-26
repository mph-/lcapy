from .expr import equation
from .nexpr import DiscreteTimeDomainExpression, nexpr

class DifferenceEquation(DiscreteTimeDomainExpression):

    def __init__(self, lhs, rhs, inputsym='x', outputsym='y', **assumptions):

        lhs = nexpr(lhs, **assumptions)
        rhs = nexpr(rhs, **assumptions)        
        
        lhssymbols = lhs.symbols
        rhssymbols = rhs.symbols

        if inputsym not in rhssymbols:
            raise ValueError('Input symbol %s not in rhs %s' % (inputsym, rhs))
        if outputsym not in lhssymbols:
            raise ValueError('Output symbol %s not in lhs %s' % (outputsym, lhs))

        super (DifferenceEquation, self).__init__(equation(lhs, rhs, **assumptions))
        self.inputsym = inputsym
        self.outputsym = outputsym        

    def transfer_function(self):
        """Create discrete-time transfer function."""
        from .zexpr import zexpr
        
        X = self.inputsym.upper()        
        Y = self.outputsym.upper()

        result = self.ZT().solve(Y)[0] / zexpr(X + '(z)')
        return result
    
    @property
    def lhs(self):
        return DiscreteTimeDomainExpression(self.expr.lhs, **self.assumptions)

    @property
    def rhs(self):
        return DiscreteTimeDomainExpression(self.expr.rhs, **self.assumptions)
    
    def dti_filter(self):
        """Create linear discrete-time invariant filter."""

        return self.transfer_function().dti_filter()
    


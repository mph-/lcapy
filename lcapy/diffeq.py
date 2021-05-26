from .expr import equation
from .nexpr import DiscreteTimeDomainExpression, nexpr

class DifferenceEquation(DiscreteTimeDomainExpression):

    def __init__(self, lhs, rhs, inputsym='x', outputsym='y', **assumptions):

        lhssymbols = lhs.symbols
        rhssymbols = rhs.symbols

        if inputsym not in rhssymbols:
            raise ValueError('Input symbol %s not in rhs %s' % (inputsym, rhs))
        if outputsym not in lhssymbols:
            raise ValueError('Output symbol %s not in lhs %s' % (outputsym, lhs))
        if inputsym in lhssymbols:
            raise ValueError('Input symbol %s in lhs %s' % (inputsym, lhs))
        if outputsym in rhssymbols:
            raise ValueError('Output symbol %s in rhs %s' % (outputsym, rhs))

        super (DifferenceEquation, self).__init__(equation(lhs, rhs, **assumptions))
        self.inputsym = inputsym
        self.outputsym = outputsym        

    def transfer_function(self):

        x = nexpr('%s(n)' % self.inputsym)
        y = nexpr('%s(n)' % self.outputsym)

        X = x.ZT()
        Y = y.ZT()

        # FIXME, it may be better to solve rather than rely on simplify
        return (self.rhs.ZT() * Y / (self.lhs.ZT() * X)).simplify()
    
    @property
    def lhs(self):
        return DiscreteTimeDomainExpression(self.expr.lhs, **self.assumptions)

    @property
    def rhs(self):
        return DiscreteTimeDomainExpression(self.expr.rhs, **self.assumptions)
    
    
    
    

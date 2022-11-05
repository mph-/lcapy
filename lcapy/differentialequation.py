"""This module provides differential equation support.

Copyright 2021 Michael Hayes, UCECE

"""

from .expr import equation, expr
from .texpr import TimeDomainExpression, texpr
import sympy as sym


class DifferentialEquation(TimeDomainExpression):

    def __init__(self, lhs, rhs, inputsym='x', outputsym='y', **assumptions):

        lhs = texpr(lhs, **assumptions)
        rhs = texpr(rhs, **assumptions)

        lhssymbols = lhs.symbols
        rhssymbols = rhs.symbols

        if inputsym not in rhssymbols:
            raise ValueError('Input symbol %s not in rhs %s' % (inputsym, rhs))
        if outputsym not in lhssymbols:
            raise ValueError('Output symbol %s not in lhs %s' %
                             (outputsym, lhs))

        cls = lhs.__class__
        eqn = sym.Eq(lhs.expr, rhs.expr, evaluate=False)
        super(DifferentialEquation, self).__init__(cls(eqn), **assumptions)
        self.inputsym = inputsym
        self.outputsym = outputsym

    def __eq__(self, x):

        if not isinstance(x, self.__class__):
            return False

        return self.lhs == x.lhs and self.rhs == x.rhs

    def transfer_function(self):
        """Create continuous-time transfer function assuming zero initial
        conditions."""
        from .sexpr import sexpr

        X = self.inputsym.upper()
        Y = self.outputsym.upper()

        Ys = self.LT(zero_initial_conditions=True).solve(Y)[0]
        Xs = sexpr(X + '(s)')

        result = Ys / Xs
        result.is_causal = True

        return result

    @property
    def lhs(self):
        return TimeDomainExpression(self.expr.lhs, **self.assumptions)

    @property
    def rhs(self):
        return TimeDomainExpression(self.expr.rhs, **self.assumptions)

    def lti_filter(self):
        """Create continuous-time linear time-invariant filter."""

        return self.transfer_function().lti_filter()

    def separate(self):
        """Rewrite differential equation so that input symbols are on the right
        and output symbols are on the left."""

        newlhs = 0
        newrhs = 0

        for term in self.lhs.as_ordered_terms():
            symbols = texpr(term).symbols
            if self.outputsym in symbols:
                newlhs += term
            if self.inputsym in symbols:
                newrhs -= term

        for term in self.rhs.as_ordered_terms():
            symbols = texpr(term).symbols
            if self.outputsym in symbols:
                newlhs -= term
            if self.inputsym in symbols:
                newrhs += term

        return self.__class__(newlhs, newrhs, **self.assumptions)

"""This module provides difference equation support.

Copyright 2021 Michael Hayes, UCECE

"""

from .expr import equation, expr
from .nexpr import DiscreteTimeDomainExpression, nexpr
import sympy as sym


class DifferenceEquation(DiscreteTimeDomainExpression):

    def __init__(self, lhs, rhs, inputsym='x', outputsym='y', **assumptions):

        lhs = nexpr(lhs, **assumptions)
        rhs = nexpr(rhs, **assumptions)

        lhssymbols = lhs.symbols
        rhssymbols = rhs.symbols

        if inputsym not in rhssymbols:
            raise ValueError('Input symbol %s not in rhs %s' % (inputsym, rhs))
        if outputsym not in lhssymbols:
            raise ValueError('Output symbol %s not in lhs %s' %
                             (outputsym, lhs))

        cls = lhs.__class__
        eqn = sym.Eq(lhs.expr, rhs.expr, evaluate=False)
        super(DifferenceEquation, self).__init__(cls(eqn), **assumptions)
        self.inputsym = inputsym
        self.outputsym = outputsym

    def __eq__(self, x):

        if not isinstance(x, self.__class__):
            return False

        return self.lhs == x.lhs and self.rhs == x.rhs

    def transfer_function(self):
        """Create discrete-time transfer function."""
        from .zexpr import zexpr

        X = self.inputsym.upper()
        Y = self.outputsym.upper()

        result = self.ZT().solve(Y)[0] / zexpr(X + '(z)')
        result.is_causal = True

        return result

    @property
    def lhs(self):
        return DiscreteTimeDomainExpression(self.expr.lhs, **self.assumptions)

    @property
    def rhs(self):
        return DiscreteTimeDomainExpression(self.expr.rhs, **self.assumptions)

    def dlti_filter(self):
        """Create discrete-time linear time-invariant filter."""

        return self.transfer_function().dlti_filter()

    def separate(self):
        """Rewrite difference equation so that input symbols are on the right
        and output symbols are on the left."""

        newlhs = 0
        newrhs = 0

        for term in self.lhs.as_ordered_terms():
            symbols = nexpr(term).symbols
            if self.outputsym in symbols:
                newlhs += term
            if self.inputsym in symbols:
                newrhs -= term

        for term in self.rhs.as_ordered_terms():
            symbols = nexpr(term).symbols
            if self.outputsym in symbols:
                newlhs -= term
            if self.inputsym in symbols:
                newrhs += term

        return self.__class__(newlhs, newrhs, **self.assumptions)


def difference_equation(lhs, rhs, inputsym='x', outputsym='y', **assumptions):
    """Create an Lcapy difference equation.

    This is an Lcapy expression of the form Eq(lhs, rhs).
    For example,
    e = difference_equation('y(n)', 'x(n) + 2 * y(n - 1)')

    The left hand side (lhs) and right hand side subexpressions
    can be obtained with the `lhs` and `rhs` attributes."""

    from .differenceequation import DifferenceEquation

    lhs = expr(lhs)
    rhs = expr(rhs)
    # Check if lhs and rhs compatible.
    diff = lhs - rhs

    return DifferenceEquation(lhs, rhs, inputsym, outputsym, **assumptions)

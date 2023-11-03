"""This module provides the Equation class.

Copyright 2023 Michael Hayes, UCECE

"""

from .expr import ExprPrint, ExprMisc, expr
from .printing import latex
from sympy import Eq


class Equation(ExprPrint, ExprMisc):

    is_Equality = True

    def __init__(self, lhs, rhs):

        self.lhs = expr(lhs)
        self.rhs = expr(rhs)

    @property
    def _pexpr(self):
        """Return expression for printing."""

        return self

    # TODO, perhaps add canonical, etc.

    def __call__(self, arg, **assumptions):
        """Transform domain."""

        from .transform import call

        return Equation(call(self.lhs, arg, **assumptions),
                        call(self.rhs, arg, **assumptions))

    @property
    def symbols(self):
        """Return dictionary of symbols in the equation keyed by name."""

        return dict(self.lhs.symbols, **self.rhs.symbols)

    def _repr_latex_(self):
        """This is used by jupyter notebooks to display an expression using
        LaTeX markup.  However, this requires mathjax.  If this method
        is not defined, jupyter falls back on _repr_pretty_ which
        outputs unicode."""

        # This is called for Expr but not ExprList
        s = latex(self._pexpr, mode='plain')
        return "$$%s$$" % s

"""This module provides the Equation class.

Copyright 2023 Michael Hayes, UCECE

"""

from .expr import ExprPrint, ExprMisc, expr
from sympy import Eq


class Equation(ExprPrint, ExprMisc):

    is_Equality = True

    def __init__(self, lhs, rhs):

        self.lhs = expr(lhs)
        self.rhs = expr(rhs)

    @property
    def _pexpr(self):
        """Return expression for printing."""

        return Eq(self.lhs.sympy, self.rhs.sympy)

    # TODO add transform

    # TODO, perhaps add canonical, etc.

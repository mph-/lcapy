"""This module provides the UndefinedDomainExpression class.

Copyright 2022 Michael Hayes, UCECE

"""

from .expr import Expr, symbol


class UndefinedDomainExpression(Expr):

    def __init__(self, val, var, **assumptions):

        self.var = symbol(var).sympy
        assumptions['var'] = self.var

        super(UndefinedDomainExpression, self).__init__(val, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

expressionclasses.register('undefined', UndefinedDomainExpression)

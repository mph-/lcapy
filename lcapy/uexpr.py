"""This module provides the UndefinedDomainExpression class.

Copyright 2022 Michael Hayes, UCECE

"""

from .expr import Expr, expr_make
from .domains import UndefinedDomain
from .sym import usersymbol
from warnings import warn

__all__ = ('uexpr', )


class UndefinedDomainExpression(UndefinedDomain, Expr):

    def __init__(self, val, var, **assumptions):

        self.var = usersymbol(var, override=False)
        assumptions['var'] = self.var

        super(UndefinedDomainExpression, self).__init__(val, **assumptions)

        if not self.has(self.var):
            warn('Expression %s does not contain variable %s' % (self, self.var))


def uexpr(arg, var, **assumptions):
    """Create UndefinedDomainExpression object."""

    assumptions['var'] = var
    return expr_make('undefined', arg, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

expressionclasses.register('undefined', UndefinedDomainExpression)

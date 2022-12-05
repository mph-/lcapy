"""This module provides the UndefinedDomainExpression class.

Copyright 2022 Michael Hayes, UCECE

"""

from .expr import Expr, expr_make
from .domains import UndefinedDomain
from .sym import usersymbol
from warnings import warn

__all__ = ('uexpr', )


class UndefinedDomainExpression(UndefinedDomain, Expr):

    def __init__(self, val, var=None, **assumptions):

        if var is None and isinstance(val, Expr):
            var = val.var
        if var is None:
            raise ValueError('Unspecified independent variable')

        self.var = usersymbol(var, override=False)
        assumptions['var'] = self.var

        super(UndefinedDomainExpression, self).__init__(val, **assumptions)

        if not self.has(self.var):
            warn('Expression %s does not contain variable %s' % (self, self.var))

    def plot(self, x=None, xlabel='', **kwargs):
        """Plot the expression.  If x is not specified, it defaults to the
        range(-0.2, 2).  x can be a vector of specified instants, a
        tuple specifing the range, or a constant specifying the
        maximum value with the minimum value set to 0.

        kwargs include:
        `axes` - the plot axes to use otherwise a new figure is created
        `xlabel` - the x-axis label
        `ylabel` - the y-axis label
        `xscale` - the x-axis scaling, say for plotting as ms
        `yscale` - the y-axis scaling, say for plotting mV
        in addition to those supported by the matplotlib plot command.

        The plot axes are returned."""

        from .plot import plot_time
        return plot_time(self, x, xlabel=xlabel, **kwargs)


def uexpr(arg, var, **assumptions):
    """Create UndefinedDomainExpression object where `var` is the
    independent variable."""

    assumptions['var'] = var
    return expr_make('undefined', arg, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

expressionclasses.register('undefined', UndefinedDomainExpression)

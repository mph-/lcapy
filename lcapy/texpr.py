"""This module provides the TimeDomainExpression class to represent
time domain expressions.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import TimeDomain
from .expr import Expr
from .functions import exp
from .sym import fsym, ssym, tsym, j, oo, tausym
from .laplace import laplace_transform
from .fourier import fourier_transform
from .units import u as uu
from sympy import Heaviside, limit, Integral, Expr as symExpr


__all__ = ('texpr', )

class TimeDomainExpression(TimeDomain, Expr):
    """t-domain expression or symbol."""

    var = tsym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)        
        assumptions['real'] = True
        super(TimeDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr        
        if check and expr.find(ssym) != set() and not expr.has(Integral):
            raise ValueError(
                't-domain expression %s cannot depend on s' % expr)
        if check and expr.find(fsym) != set() and not expr.has(Integral):
            raise ValueError(
                't-domain expression %s cannot depend on f' % expr)

    def _mul_compatible_domains(self, x):

        if self.domain == x.domain:
            return True        

        return x.is_constant_domain

    def _div_compatible_domains(self, x):

        if self.domain == x.domain:
            return True
        
        return x.is_constant_domain

    def as_expr(self):
        return TimeDomainExpression(self)

    def infer_assumptions(self):

        self.assumptions.infer_from_expr(self)

    def LT(self, evaluate=True, **assumptions):
        """Determine one-sided Laplace transform with 0- as the lower limit.

        This is an alias for laplace."""

        return self.laplace(evaluate, **assumptions)
            
    def laplace(self, evaluate=True, **assumptions):
        """Determine one-sided Laplace transform with 0- as the lower limit."""

        assumptions = self.assumptions.merge_and_infer(self, **assumptions)
        result = laplace_transform(self.expr, self.var, ssym, evaluate=evaluate)
        return self.change(result, domain='laplace', units_scale=uu.s, **assumptions)

    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        from .phasor import PhasorTimeDomainExpression
        
        return PhasorTimeDomainExpression.from_time(self, **assumptions)

    def FT(self, evaluate=True, **assumptions):
        """Attempt Fourier transform.  This is an alias for fourier.

        X(f) = \int_{-\infty}^{\infty} x(t) exp(-j 2\pi f t) dt."""

        return self.fourier(evaluate, **assumptions)        

    def fourier(self, evaluate=True, **assumptions):
        """Attempt Fourier transform."""

        assumptions = self.assumptions.merge_and_infer(self, **assumptions)
        result = fourier_transform(self.expr, self.var, fsym, evaluate=evaluate)
        return self.change(result, domain='fourier', units_scale=uu.s, **assumptions)

    def angular_fourier(self, evaluate=True, **assumptions):
        """Attempt angular Fourier transform."""

        from .symbols import omega, pi, f

        assumptions = self.assumptions.merge_and_infer(self, **assumptions)
        result = self.fourier(evaluate, **assumptions).subs(f, omega / (2 * pi))
        # Could optimise...
        return self.change(result, domain='angular fourier', units_scale=uu.s, **assumptions)

    def time(self, **assumptions):
        return self

    def plot(self, t=None, **kwargs):
        """Plot the time waveform.  If t is not specified, it defaults to the
        range (-0.2, 2).  t can be a vector of specified instants, a
        tuple specifing the range, or a constant specifying the
        maximum value with the minimum value set to 0.

        kwargs include:
        axes - the plot axes to use otherwise a new figure is created
        xlabel - the x-axis label
        ylabel - the y-axis label
        xscale - the x-axis scaling, say for plotting as ms
        yscale - the y-axis scaling, say for plotting mV
        in addition to those supported by the matplotlib plot command.
        
        The plot axes are returned."""

        from .plot import plot_time
        return plot_time(self, t, **kwargs)

    def sample(self, t):

        """Return a discrete-time signal evaluated at time values specified by
        vector t. """

        return self.evaluate(t)

    def initial_value(self):
        """Determine value at t = 0. 
        See also pre_initial_value and post_initial_value"""

        return self.subs(0)

    def pre_initial_value(self):
        """Determine value at t = 0-.
        See also initial_value and post_initial_value"""

        return self.limit(self.var, 0, dir='-')

    def post_initial_value(self):
        """Determine value at t = 0+.
        See also pre_initial_value and initial_value"""

        return self.limit(self.var, 0, dir='+')    

    def final_value(self):
        """Determine value at t = oo."""

        return self.limit(self.var, oo)

    def remove_condition(self):
        """Remove the piecewise condition from the expression.
        See also force_causal."""

        if not self.is_conditional:
            return self
        expr = self.expr
        expr = expr.args[0].args[0]
        return self.__class__(expr)

    def force_causal(self):
        """Remove the piecewise condition from the expression
        and multiply by Heaviside function.  See also remove_condition."""

        if self.is_causal:
            return self
        
        expr = self.expr
        if self.is_conditional:
            expr = expr.args[0].args[0]            
        expr = expr * Heaviside(t)
        return self.__class__(expr)        

    def convolve(self, impulseresponse, commutate=False, **assumptions):
        """Convolve self with impulse response."""

        if not isinstance(impulseresponse, TimeDomainExpression):
            raise ValueError('Expecting TimeDomainExpression for impulse response')

        f1 = self.expr
        f2 = impulseresponse.expr
        if commutate:
            f1, f2 = f2, f1
        result = Integral(f1.subs(self.var, self.var - tausym) *
                          f2.subs(self.var, tausym),
                          (tausym, -oo, oo))
        return self.__class__(result, **assumptions)

    
class TimeDomainImpulseResponse(TimeDomainExpression):
    """Time-domain impulse response."""

    # TODO, check attributes.
    quantity = 'transfer'
    quantity_label = 'Impulse response'
    domain_units = '1/s'
    is_transfer = True


def texpr(arg, **assumptions):
    """Create TimeDomainExpression object.  If `arg` is tsym return t"""

    if arg is tsym:
        return t
    return TimeDomainExpression(arg, **assumptions)


from .expressionclasses import expressionclasses

classes = expressionclasses.register('time', TimeDomainExpression)
TimeDomainVoltage = classes['voltage']
TimeDomainCurrent = classes['current']

t = TimeDomainExpression('t')
t.units = uu.s

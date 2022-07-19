"""This module provides the TimeDomainExpression class to represent
time domain expressions.

Copyright 2014--2022 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import TimeDomain
from .expr import Expr, expr_make
from .functions import exp
from .state import state
from .sym import fsym, ssym, tsym, j, oo
from .laplace import laplace_transform
from .fourier import fourier_transform
from .units import u as uu
from sympy import Heaviside, Derivative, Integral, limit, Expr as symExpr


__all__ = ('texpr', )


class TimeDomainExpression(TimeDomain, Expr):
    """t-domain expression or symbol."""

    var = tsym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        assumptions['real'] = True
        super(TimeDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr
        # Should be more specific...
        if check and expr.has(ssym) and not expr.has(Integral) and not expr.has(Derivative):
            raise ValueError(
                't-domain expression %s cannot depend on s' % expr)
        if check and expr.has(fsym) and not expr.has(Integral):
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

    @property
    def abs(self):
        """Return absolute value."""

        return self.__class__(abs(self.expr), **self.assumptions)

    def as_expr(self):
        return TimeDomainExpression(self)

    def infer_assumptions(self):

        self.assumptions.infer_from_expr(self)

    def LT(self, evaluate=True, zero_initial_conditions=None, **assumptions):
        """Determine one-sided Laplace transform with 0- as the lower limit.

        By default initial conditions are assumed to be zero.  This
        can be controlled by `zero_initial_conditions`.

        This is an alias for laplace.

        """
        if zero_initial_conditions is None:
            zero_initial_conditions = state.zero_initial_conditions

        assumptions = self.assumptions.merge_and_infer(self, **assumptions)
        result = laplace_transform(
            self.expr, self.var, ssym, evaluate=evaluate,
            zero_initial_conditions=zero_initial_conditions)
        return self.change(result, domain='laplace', units_scale=uu.s, **assumptions)

    def laplace(self, evaluate=True, zero_initial_conditions=None,
                **assumptions):
        """Determine one-sided Laplace transform with 0- as the lower limit.

        By default initial conditions are assumed to be zero.  This
        can be controlled by `zero_initial_conditions`.

        This is an alias for LT."""

        return self.LT(evaluate=evaluate, zero_initial_conditions=zero_initial_conditions,
                       **assumptions)

    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        from .phasor import PhasorTimeDomainExpression

        return PhasorTimeDomainExpression.from_time(self, **assumptions)

    def FT(self, var=None, evaluate=True, **assumptions):
        """Attempt Fourier transform.

        X(f) = \int_{-\infty} ^ {\infty} x(t) exp(-j 2\pi f t) dt."""

        from .symbols import f, omega, Omega, F

        if var is None:
            var = f
        if id(var) not in (id(f), id(F), id(omega), id(Omega)):
            raise ValueError(
                'FT requires var to be f, F, omega, or Omega`, not %s' % var)

        assumptions = self.assumptions.merge_and_infer(self, **assumptions)
        result = fourier_transform(
            self.expr, self.var, fsym, evaluate=evaluate)
        result = self.change(result, domain='fourier',
                             units_scale=uu.s, **assumptions)
        result = result(var)
        result = result.expand(diracdelta=True, wrt=var)
        result = result.simplify()
        return result

    def fourier(self, var=None, evaluate=True, **assumptions):
        """Attempt Fourier transform. This is an alias for FT."""

        return self.FT(var, evaluate, **assumptions)

    def norm_fourier(self, evaluate=True, **assumptions):
        """Attempt normalized Fourier transform."""

        from .symbols import F

        return self.FT(F, evaluate, **assumptions)

    def angular_fourier(self, evaluate=True, **assumptions):
        """Attempt angular Fourier transform."""

        from .symbols import omega

        if self.has(omega):
            raise ValueError(
                'Time-domain expression has omega.  Substitute with another symbol, say omega0')

        return self.FT(omega, evaluate, **assumptions)

    def norm_angular_fourier(self, evaluate=True, **assumptions):
        """Attempt normalized angular Fourier transform."""

        from .symbols import Omega

        return self.FT(Omega, evaluate, **assumptions)

    def time(self, **assumptions):
        return self

    def plot(self, t=None, **kwargs):
        """Plot the time waveform.  If t is not specified, it defaults to the
        range(-0.2, 2).  t can be a vector of specified instants, a
        tuple specifing the range, or a constant specifying the
        maximum value with the minimum value set to 0.

        kwargs include:
        `axes` - the plot axes to use otherwise a new figure is created
        `xlabel` - the x-axis label
        `ylabel` - the y-axis label
        `xscale` - the x-axis scaling, say for plotting as ms
        `yscale` - the y-axis scaling, say for plotting mV
        `plot_deltas` - plot Dirac deltas as arrows
        in addition to those supported by the matplotlib plot command.

        The plot axes are returned."""

        from .plot import plot_time
        return plot_time(self, t, **kwargs)

    def response(self, xvector, tvector, method='bilinear'):
        """Evaluate response to input signal `xvector` at times
        `tvector`.  This returns a NumPy array."""

        from .sexpr import s
        return self(s).response(xvector, tvector, method=method)

    def sample(self, tvector):
        """Return a discrete-time signal evaluated at time values
        specified by `tvector`.  This returns a NumPy array."""

        return self.evaluate(tvector)

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

    def force_causal(self):
        """Remove the piecewise condition from the expression
        and multiply by Heaviside function.  See also remove_condition."""

        if self.is_causal:
            return self

        expr = self.expr
        if self.is_conditional:
            expr = expr.args[0].args[0]
        expr = expr * Heaviside(t.var)
        return self.__class__(expr)

    def discretize(self, method='bilinear', alpha=0.5):
        """Convert to a discrete-time approximation:

        : math: `h[n] \approx h(t)`

        With the impulse-invariance method, the discrete-time impulse response
        is related to the continuous-time impulse response by

        : math: `h[n] = h_c(n \Delta t)`

        Note, when designing digital filters, it is often common to to
        scale the discrete-time impulse response by the sampling
        interval:

        : math: `h[n] = \Delta t h_c(n \Delta t)`

        The default method is 'bilinear'.  Other methods are:
        'impulse-invariance' 'bilinear', 'tustin', 'trapezoidal'
        'generalized-bilinear', 'gbf' controlled by the parameter
        `alpha` 'euler', 'forward-diff', 'forward-euler'
        'backward-diff', 'backward-euler' 'simpson', 'matched-Z',
        'zero-pole-matching'

        """

        from .sym import dt
        from .symbols import n

        if method in ('impulse-invariance', ):
            return self.subs(t, n * dt)

        Hc = self.LT()
        H = Hc.discretize(method=method, alpha=alpha)
        return H.IZT()

    def zdomain(self, **assumptions):

        # Going via the Laplace domain will restrict result for n >= 0.
        # Perhaps just sample with t = n * dt and warn if expression
        # contains Dirac deltas?
        return self.laplace().discretize(**assumptions)

    def discrete_frequency(self, **assumptions):

        return self.discrete_time(**assumptions).discrete_frequency()

    def discrete_time(self, **assumptions):

        return self.discretize(**assumptions)

    def differential_equation(self, inputsym='x', outputsym='y'):
        """Create differential equation from impulse response.
        """

        H = self.LT(zero_initial_conditions=True)
        return H.differential_equation(inputsym, outputsym)


class TimeDomainImpulseResponse(TimeDomainExpression):
    """Time-domain impulse response."""

    # TODO, check attributes.
    quantity = 'transfer'
    quantity_label = 'Impulse response'
    domain_units = uu.Hz
    is_transfer = True


def texpr(arg, **assumptions):
    """Create TimeDomainExpression object.  If `arg` is tsym return t"""

    if arg is tsym:
        return t
    return expr_make('time', arg, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

classes = expressionclasses.register('time', TimeDomainExpression)
TimeDomainVoltage = classes['voltage']
TimeDomainCurrent = classes['current']

t = TimeDomainExpression('t')
t.units = uu.s

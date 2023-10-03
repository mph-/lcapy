"""This module provides the TimeDomainExpression class to represent
time domain expressions.

Copyright 2014--2023 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import TimeDomain
from .expr import Expr, expr_make
from .functions import exp, integrate
from .state import state, validate
from .sym import fsym, omegasym, ssym, tsym, j, oo
from .laplace import laplace_transform
from .fourier import fourier_transform
from .units import u as uu
from sympy import Heaviside, limit, Expr as symExpr


__all__ = ('TimeDomainExpression', 'texpr', )


class TimeDomainExpression(TimeDomain, Expr):
    """t-domain expression or symbol."""

    var = tsym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        assumptions['real'] = True
        super(TimeDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr

        if check and not expr.has(tsym):
            if expr.has(fsym):
                validate(state.f_in_t,
                         't-domain expression %s depends on f' % expr)
            if expr.has(ssym):
                validate(state.s_in_t,
                         't-domain expression %s depends on s' % expr)
            if expr.has(omegasym):
                validate(state.w_in_t,
                         't-domain expression %s depends on omega' % expr)

    def _mul_compatible_domains(self, x):

        if self.domain == x.domain:
            return True

        return x.is_constant_domain

    def _mul_dubious(self, x):

        # Allow rect(t) * rect(t) etc.
        if not (self.is_signal or self.is_immittance or
                x.is_signal or x.is_immittance):
            return False

        # Allow v*v.
        return x.is_time_domain and not (self.is_signal and x.is_signal)

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

    @property
    def energy(self):
        """Return signal energy."""

        # TODO: teach SymPy improper integral of rect and sinc functions

        # Only need abs for complex signal.
        E = integrate(abs(self)**2, (t, -oo, oo))

        # FIXME: result should not be a time-domain expression
        return E

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

        if self.has(ssym):
            raise ValueError(
                'Time-domain expression has s.  Substitute with another symbol, say p')

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

        from .phasor import PhasorDomainExpression

        return PhasorDomainExpression.from_time(self, **assumptions)

    def FT(self, var=None, evaluate=True, **assumptions):
        """Attempt Fourier transform.

        X(f) = \int_{-\infty} ^ {\infty} x(t) exp(-j 2\pi f t) dt."""

        from .symbols import f, omega, Omega, F

        if var is None:
            var = f
        if id(var) not in (id(f), id(F), id(omega), id(Omega)):
            raise ValueError(
                'FT requires var to be f, F, omega, or Omega`, not %s' % var)

        if self.has(var):
            raise ValueError(
                'Time-domain expression has %s.  Substitute with another symbol, say %s0' % (var, var))

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

        return self.FT(omega, evaluate, **assumptions)

    def norm_angular_fourier(self, evaluate=True, **assumptions):
        """Attempt normalized angular Fourier transform."""

        from .symbols import Omega

        return self.FT(Omega, evaluate, **assumptions)

    def frequency_response(self, **assumptions):
        """Convert to frequency response domain.  Note, this is similar to the
        Fourier domain but not always."""

        return self.laplace(**assumptions).frequency_response()

    def angular_frequency_response(self, **assumptions):
        """Convert to angular frequency response domain.  Note, this is
        similar to the angular Fourier domain but not always."""

        return self.laplace(**assumptions).angular_frequency_response()

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

    def discretize(self, method=None, alpha=0.5):
        """Convert to a discrete-time approximation:

        :math:`x(n) \approx K x(t)`

        where :math:`K` is a scale factor.

        The default method depends on the expression quantity.  If the
        quantity is `undefined`, `voltage`, `current,
        `voltagesquared`, `currentsquared`, or `power`, the default
        method is `impulse-invariance` and there is no scaling.  This
        uses:

        : math: `x[n] = x(n \Delta t)`

        For other quantities, the default method is `bilinear', and
        the result is scaled by the sampling interval (`Delta t`) as
        is common for the discrete-time impulse response of digital
        filters:

        : math: `h[n] = \Delta t h_c(n \Delta t)`

        The methods are: 'impulse-invariance' 'bilinear', 'tustin',
        'trapezoidal' 'generalized-bilinear', 'gbf' controlled by the
        parameter `alpha` 'euler', 'forward-diff', 'forward-euler'
        'backward-diff', 'backward-euler' 'simpson', 'matched-Z',
        'zero-pole-matching'.  All the methods except
        `impulse-invariance` are applied in the Laplace domain so the
        result is unknown for `n < 0` unless the expression is causal.

        """

        from .sym import dt
        from .symbols import n

        signal = self.is_signal or self.is_squared or \
            self.is_power or self.is_undefined

        if method is None:
            if signal:
                method = 'impulse-invariance'
            else:
                method = 'bilinear'

        if method in ('impulse-invariance', ):
            from .symbols import t

            result = self.subs(t, n * dt)
            if not signal:
                result *= dt
            return result

        Hc = self.LT()
        H = Hc.discretize(method=method, alpha=alpha)
        if self.is_undefined:
            H /= dt
        return H.IZT()

    def dlti_filter(self, method='bilinear'):
        """Create DLTI filter using bilinear transform."""

        return self.LT().dlti_filter(method)

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


t = TimeDomainExpression('t')
t.units = uu.s

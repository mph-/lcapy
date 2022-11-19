"""This module provides the LaplaceDomainExpression class to represent
s-domain (Laplace domain) expressions.

Copyright 2014--2022 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import LaplaceDomain
from .inverse_laplace import inverse_laplace_transform
from .state import state, validate
from .sym import ssym, tsym, fsym, omegasym, j, pi, sympify, dt
from .ratfun import _zp2tf, _pr2tf, Ratfun
from .expr import Expr, symbol, expr, ExprDict, ExprList, exprcontainer, expr_make
from .differentialequation import DifferentialEquation
from .units import u as uu
from .functions import sqrt, DiracDelta, Function
from sympy import limit, exp, Poly, Derivative, Integral, div, oo, Eq, Expr as symExpr
from warnings import warn


__all__ = ('sexpr', 'zp2tf', 'tf', 'pr2tf')


class LaplaceDomainExpression(LaplaceDomain, Expr):
    """s-domain expression or symbol."""

    var = ssym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        super(LaplaceDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr

        if check and not expr.has(ssym):
            if expr.has(fsym):
                validate(state.f_in_s,
                         's-domain expression %s depends on f' % expr)
            if expr.has(tsym):
                validate(state.t_in_s,
                         's-domain expression % s depends on t' % expr)
            if expr.has(omegasym):
                validate(state.w_in_s,
                         's-domain expression %s depends on omega' % expr)

    def as_expr(self):
        return LaplaceDomainExpression(self)

    @classmethod
    def from_poles_residues(cls, poles, residues):
        """Create a transfer function from lists of poles and residues.

        See also `from_zeros_poles_gain`, `from_numer_denom`."""

        return cls(pr2tf(poles, residues, cls.var), causal=True)

    @classmethod
    def from_zeros_poles_gain(cls, zeros, poles, K=1):
        """Create a transfer function from lists of zeros and poles,
        and from a constant gain.

        See also `from_poles_residues`, `from_numer_denom`."""

        return cls(zp2tf(zeros, poles, K, cls.var), causal=True)

    @classmethod
    def from_numer_denom(cls, numer, denom):
        """Create a transfer function from lists of the coefficient
        for the numerator and denominator.

        See also `from_zeros_poles_gain`, `from_poles_residues`."""

        return cls(tf(numer, denom, cls.var), causal=True)

    def tdifferentiate(self):
        """Differentiate in t-domain (multiply by s)."""

        return self.__class__(self.expr * self.var, **self.assumptions)

    def tintegrate(self):
        """Integrate in t-domain (divide by s)."""

        return self.__class__(self.expr / self.var, **self.assumptions)

    def delay(self, T):
        """Apply delay of T seconds by multiplying by exp(-s T)."""

        T = self.__class__(T)
        return self.__class__(self.expr * exp(-s * T))

    @property
    def jomega(self):
        """Return expression with s = j omega."""

        from .symbols import jomega
        return self.subs(self.var, jomega)

    def post_initial_value(self):
        """Determine post-initial value at t = 0+."""

        return self.__class__(limit(self.expr * self.var, self.var, oo))

    def final_value(self):
        """Determine value at t = oo."""

        return self.__class__(limit(self.expr * self.var, self.var, 0))

    def inverse_laplace(self, zero_initial_conditions=None, **assumptions):
        """Attempt inverse Laplace transform.

        If causal=True the response is zero for t < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for t < 0.
        Otherwise the result is only known for t >= 0.

        By default initial conditions are assumed to be zero.  This
        can be controlled by `zero_initial_conditions`.

        """

        return self.ILT(zero_initial_conditions=zero_initial_conditions, **assumptions)

    def ILT(self, zero_initial_conditions=None, **assumptions):
        """Attempt inverse Laplace transform.

        If causal=True the response is zero for t < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for t < 0.
        Otherwise the result is only known for t >= 0.

        By default initial conditions are assumed to be zero.  This
        can be controlled by `zero_initial_conditions`.

        """
        if zero_initial_conditions is None:
            zero_initial_conditions = state.zero_initial_conditions

        assumptions = self.assumptions.merge(**assumptions)
        result = inverse_laplace_transform(self.expr, self.var, tsym,
                                           zero_initial_conditions=zero_initial_conditions,
                                           **assumptions)

        return self.change(result, domain='time', units_scale=uu.Hz,
                           **assumptions)

    def time(self, **assumptions):
        """Convert to time domain.

        If causal=True the response is zero for t < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for t < 0.
        Otherwise the result is only known for t >= 0.

        """

        try:
            return self.inverse_laplace(**assumptions)
        except ValueError:
            return self.as_sum().inverse_laplace(**assumptions)

    def laplace(self, **assumptions):
        """Convert to s-domain."""

        assumptions = self.assumptions.merge(**assumptions)
        return self.__class__(self, **assumptions)

    def fourier(self, **assumptions):
        """Convert to Fourier domain."""
        from .symbols import f, jw, pi

        # Require f = 0 in the region of convergence.  This is
        # achieved if the expression is the Laplace transform of a
        # stable impulse response.
        if ((self.is_causal or assumptions.get('causal', False))
                and self.is_stable):
            tmp = self.subs(jw)
            return self.change(tmp.subs(2 * pi * f, safe=True),
                               domain='fourier', **assumptions)

        result = self.time(**assumptions).fourier(**assumptions)
        return result

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import jw

        if ((self.is_causal or assumptions.get('causal', False))
                and self.is_stable):
            tmp = self(jw)
            return self.change(tmp, domain='angular fourier',
                               **assumptions)

        result = self.time(**assumptions).angular_fourier(**assumptions)
        return result

    def norm_angular_fourier(self, **assumptions):
        """Convert to normalized angular Fourier domain."""
        from .symbols import jw, Omega

        if ((self.is_causal or assumptions.get('causal', False))
                and self.is_stable):
            tmp = self(j * Omega / dt)
            return self.change(tmp, domain='norm angular fourier',
                               **assumptions)

        result = self.time(**assumptions).norm_angular_fourier(**assumptions)
        return result

    def norm_fourier(self, **assumptions):
        """Convert to normalized Fourier domain."""
        from .symbols import jw, F

        if ((self.is_causal or assumptions.get('causal', False))
                and self.is_stable):
            tmp = self(j * F / dt)
            return self.change(tmp, domain='norm fourier',
                               **assumptions)

        result = self.time(**assumptions).norm_fourier(**assumptions)
        return result

    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        result = PhasorDomainExpression.from_laplace(self, **assumptions)
        return result

    def phasor_ratio(self, **assumptions):
        """Convert to phasor ratio domain."""

        result = PhasorRatioDomainExpression.from_laplace(self, **assumptions)
        return result

    def transient_response(self, tvector=None):
        """Evaluate transient (impulse) response."""

        if tvector is None:
            return self.time()

        return self.time().evaluate(tvector)

    def impulse_response(self, tvector=None):
        """Evaluate transient (impulse) response."""

        return self.transient_response(tvector)

    def step_response(self, tvector=None):
        """Evaluate step response."""

        H = self.__class__(self / self.var, **self.assumptions)
        return H.transient_response(tvector)

    def frequency_response(self, **assumptions):
        """Convert to frequency response domain.  Note, this is similar to the
        Fourier domain but not always."""

        from .symbols import jf

        tmp = self.subs(jf * 2 * pi)
        return self.change(tmp, domain='frequency response', **assumptions)

    def angular_frequency_response(self, **assumptions):
        """Convert to angular frequency response domain.  Note, this is
        similar to the angular Fourier domain but not always."""
        from .symbols import jw

        tmp = self.subs(jw)
        return self.change(tmp, domain='angular frequency response',
                           **assumptions)

    def frequency_response_evaluate(self, fvector=None, var=None,
                                    **assumptions):

        from .symbols import f, jf, omega, jw

        if var in (None, f, jf):
            X = self.frequency_response()
        elif var in (omega, jw):
            X = self.angular_frequency_response()
        else:
            raise ValueError('Invalid var %s' % var)

        if fvector is None:
            return X

        return X.evaluate(fvector)

    def _response_impulse_invariance(self, xvector, tvector, dtval):
        """Evaluate response of system with applied signal using impulse
        invariance method.  This method assumes that the data is well
        over-sampled and that there are no Dirac deltas."""

        from numpy import arange, convolve, diff, hstack

        # Perform polynomial long division so expr1 = Q + M / A
        B, A, delay, undef = self._as_B_A_delay_undef()
        if undef != 1:
            raise ValueError('Have undefined expression %s' % undef)
        Q, M = div(B, A, self.var)
        expr1 = M / A

        Nt = len(tvector)

        # Evaluate impulse response.
        th = arange(Nt) * dtval
        H = LaplaceDomainExpression(expr1, **self.assumptions)
        hvector = H.transient_response(th)

        ty = tvector
        y = convolve(xvector, hvector)[0:Nt] * dtval

        if Q:
            # Handle Dirac deltas and their derivatives.
            C = expr(Q).coeffs()
            for c in reversed(C):

                y += float(c.expr) * xvector

                xvector = diff(xvector) / dtval
                xvector = hstack((xvector, 0))

        from scipy.interpolate import interp1d

        if delay != 0.0:
            # Try linear interpolation; should oversample first...
            y = interp1d(ty, y, bounds_error=False, fill_value=0)
            y = y(tvector - float(delay))

        return y

    def _response_bilinear(self, xvector, tvector, dtval, alpha=0.5):
        """Evaluate response of system with applied signal using bilinear
         method."""

        import scipy.signal as signal
        from numpy import hstack, zeros

        expr, delay = self.as_ratfun_delay()
        Ndelay = 0
        if delay != 0:
            if delay < 0:
                raise ValueError('Need time advance %s' % delay)
            Ndelay = round(float(delay / dtval), 12)
            if not Ndelay.is_integer:
                # Require fractional delay.
                Ndelay = 0
                delay = 0
                expr = self
            else:
                Ndelay = int(Ndelay)
        if delay != 0:
            xvector = xvector[:-Ndelay]
            tvector = tvector[:-Ndelay] - delay

        if expr.has(exp):
            warn('Using second order Pade approximation for exp.')
            expr = self.approximate_exp(method='pade', order=2, numer_order=1)

        Hz = expr.generalized_bilinear_transform(alpha).subs(dt, dtval)
        fil = Hz.dlti_filter()

        a = [a1.fval for a1 in fil.a]
        b = [b1.fval for b1 in fil.b]

        y = signal.lfilter(b, a, xvector)

        if Ndelay != 0:
            y = hstack((zeros(Ndelay), y))
        return y

    def response(self, xvector, tvector, method='bilinear', alpha=0.5):
        """Evaluate response to input signal `xvector` at times
        specified by `tvector`.  This returns a NumPy array.

        The default method is 'bilinear'.  Other methods are:
        'impulse-invariance'
        'bilinear', 'tustin', 'trapezoidal'
        'generalized-bilinear', 'gbf' controlled by the parameter `alpha`
        'euler', 'forward-diff', 'forward-euler'
        'backward-diff', 'backward-euler'.

        See also `discretize` with regards to scaling of the result.
        """

        from numpy import array, diff, allclose

        xvector = array(xvector)
        tvector = array(tvector)

        symbols = self.symbols
        symbols.pop('s', None)
        if symbols != {}:
            raise ValueError('Have undefined symbols: %s' % symbols)

        if len(xvector) != len(tvector):
            raise ValueError('xvector, tvector must have at least 2 samples')

        if len(xvector) < 2:
            raise ValueError('xvector must have same length as tvector')

        td = diff(tvector)
        if not allclose(diff(td), 0):
            raise ValueError('Time samples not uniformly spaced')

        dtval = td[0]
        if dtval < 0:
            raise ValueError('Time samples not increasing')

        if self.is_constant:
            return float(self.expr) * xvector

        # Expand cosh, sinh into sum of exps.
        expr = self.expand_hyperbolic_trig()

        if method in ('impulse-invariance', 'adhoc'):
            result = expr._response_impulse_invariance(xvector, tvector, dtval)
        elif method in ('bilinear', 'tustin', 'trapezoidal'):
            result = expr._response_bilinear(
                xvector, tvector, dtval, alpha=0.5)
        elif method in ('gbf', 'generalized-bilinear'):
            result = expr._response_bilinear(
                xvector, tvector, dtval, alpha=alpha)
        elif method in ('euler', 'forward-diff', 'forward-euler'):
            result = expr._response_bilinear(xvector, tvector, dtval, alpha=0)
        elif method in ('backward-diff', 'backward-euler'):
            result = expr._response_bilinear(xvector, tvector, dtval, alpha=1)
        else:
            raise ValueError('Unknown method %s' % method)

        return result * dtval

    def state_space(self, form='CCF'):
        """Create state-space representation from transfer function.  Note,
        state-space representations are not unique and are determined
        by the `form` argument.  Currently this can be 'CCF' for the
        controllable canonical form, 'OCF' for the observable
        canonical form, or 'DCF' for the diagonal canonical form."""

        from .statespace import StateSpace

        a = self.a
        b = self.b

        return StateSpace.from_transfer_function_coeffs(b, a, form)

    @property
    def ss(self):
        """Return state-space representation using controllable canonical form.
        For other forms, use `state_space()`."""

        return self.state_space()

    def differential_equation(self, inputsym='x', outputsym='y'):
        """Create differential equation from transfer function.

        For example,
        >>> H = (s + 3) / (s**2 + 4)
        >>> H.differential_equation()
                 d                    d
        3.y(t) + --(y(t)) = 4.x(t) + ---(x(t))
                 dt                    2
                                     dt
        """

        from .symbols import t

        x = Function(inputsym)(t)
        y = Function(outputsym)(t)

        X = x.LT()
        Y = y.LT()

        H = self.simplify()
        N = H.N
        D = H.D

        lhs = (N * Y).ILT(causal=True, zero_initial_conditions=True)
        rhs = (D * X).ILT(causal=True, zero_initial_conditions=True)

        return DifferentialEquation(lhs, rhs, inputsym, outputsym)

    def lti_filter(self, normalize_a0=True):
        """Create continuous-time linear time-invariant filter from
        continuous-time transfer function."""

        # TODO, perhaps add only to TimeDomainTransfer
        # or check that is_transfer, is_impedance, or is_admittance

        from .ltifilter import LTIFilter

        if not self.is_rational_function:
            raise ValueError("Not a rational function")

        N = self.N
        D = self.D
        bn = N.coeffs()
        an = D.coeffs()

        if normalize_a0:
            bn = [bx / an[0] for bx in bn]
            an = [ax / an[0] for ax in an]

        lpf = LTIFilter(bn, an)
        return lpf

    def dlti_filter(self, method='bilinear', alpha=0.5):
        """Create DLTI filter using bilinear transform."""

        from .transfer import transfer

        e = self

        if not e.is_causal:
            warn('Assuming expression is causal')
            e = expr(e, causal=True)

        if not e.is_ratio:
            warn('Assuming %s expression is a transfer function' %
                 e.quantity)
            e = transfer(e)

        return e.discretize(method, alpha).simplify().dlti_filter()

    def evaluate(self, svector=None):

        return super(LaplaceDomainExpression, self).evaluate(svector)

    def plot(self, **kwargs):
        """Plot pole-zero map.

        kwargs include:
        axes - the plot axes to use otherwise a new figure is created
        xlabel - the x-axis label (default Re(s))
        ylabel - the y-axis label (default Im(s))
        xscale - the x-axis scaling
        yscale - the y-axis scaling
        in addition to those supported by the matplotlib plot command.

        The plot axes are returned."""

        from .plot import plot_pole_zero
        return plot_pole_zero(self, **kwargs)

    def pole_zero_plot(self, **kwargs):
        """Plot pole-zero map."""

        return self.plot(**kwargs)

    def bode_plot(self, fvector=None, unwrap=True, var=None, strict=False,
                  **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Bode
        plot (but without the straight line approximations).

        If `var` is None or `f`, `fvector` specifies the (linear)
        frequencies (in hertz) otherwise if `var` is `omega`,
        `fvector` specifies the angular frequencies (in radians).

         If `fvector` is a tuple (f1, f2), it sets the frequency
        limits.  Since a logarithmic frequency scale is used, f1 must
        be greater than 0.

        `unwrap` controls phase unwrapping (default True).

        This method makes the assumption that the expression is causal.
        Note, the Bode plot does not exist for marginally stable and unstable
        systems since `jw` is outside the region of convergence.

        """

        H = self.frequency_response_evaluate(var=var)
        return H.bode_plot(fvector, unwrap=unwrap, **kwargs)

    def nyquist_plot(self, fvector=None, var=None, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Nyquist
        plot (imaginary part versus real part).  fvector specifies the
        frequencies.  If it is a tuple (f1, f2), it sets the frequency
        limits.

        `npoints` set the number of plotted points.

        The unit circle is shown by default.  This can be disabled with `unitcircle=False`.

        This method makes the assumption that the expression is causal.
        """

        H = self.frequency_response_evaluate(var=var)
        return H.nyquist_plot(fvector, **kwargs)

    def nichols_plot(self, fvector=None, var=None, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Nichols
        plot (dB versus phase).  fvector specifies the frequencies.
        If it is a tuple (f1, f2), it sets the frequency limits.

        `npoints` set the number of plotted points.

        This method makes the assumption that the expression is causal.

        """

        H = self.frequency_response_evaluate(var=var)
        return H.nichols_plot(fvector, **kwargs)

    def generalized_bilinear_transform(self, alpha=0.5):
        """Approximate s = ln(z) / dt

        by s = (1 / dt) * (1 - z**-1) / (alpha + (1 - alpha) * z**-1))

        When `alpha = 0` this is equivalent to the forward Euler
        (forward difference) method.

        When `alpha = 0.5` this is equivalent to the bilinear
        transform, aka, Tustin's method or the trapezoidal method.

        When `alpha = 1` this is equivalent to the backward Euler
        (backward difference) method.

        See also `discretize` with regards to scaling of the result.

        """

        from .discretetime import z, dt

        if alpha < 0 or alpha > 1:
            raise ValueError("alpha must be between 0 and 1 inclusive")

        scale = 1 if self.is_ratio else 1 / dt

        return self.subs((1 / dt) * (1 - z**-1) / (alpha + (1 - alpha) * z**-1)) * scale

    def bilinear_transform(self):
        """Approximate s = ln(z) / dt

        by s = (2 / dt) * (1 - z**-1) / (1 + z**-1)

        This is also called Tustin's method and is equivalent to the
        trapezoidal method.

        See also `discretize` with regards to scaling of the result.
        """

        # TODO: add frequency warping as an option
        return self.generalized_bilinear_transform(0.5)

    def forward_euler_transform(self):
        """Approximate s = ln(z) / dt

        by s = (1 / dt) * (1 - z**-1) / z**-1.   This is also known
        as the forward difference method.

        See also `discretize` with regards to scaling of the result.
        """

        return self.generalized_bilinear_transform(0)

    def backward_euler_transform(self):
        """Approximate s = ln(z) / dt

        by s = (1 / dt) * (1 - z**-1).  This is also known
        as the backward difference method.

        See also `discretize` with regards to scaling of the result.
        """

        return self.generalized_bilinear_transform(1)

    def simpson_transform(self):
        """Approximate s = ln(z) / dt

        by s = (3 / dt) * (z**2 - 1) / (z**2 + 4 * z + 1).  This is
        more accurate than the other methods but doubles the system
        order and can produce unstable poles.

        See also `discretize` with regards to scaling of the result.

        """

        from .discretetime import z

        scale = dt if self.is_ratio else 1

        return self.subs((3 / dt) * (z**2 - 1) / (z**2 + 4 * z + 1)) * scale

    def matched_ztransform(self):
        """Match poles and zeros of H(s) to approximate H(z).

        If there are no zeros, this is equivalent to impulse_invariance.

        See also `bilinear_transform` and `impulse_invariance_transform`.

        See also `discretize` with regards to scaling of the result.
        """

        from .discretetime import z

        scale = dt if self.is_ratio else 1

        zeros, poles, K, undef = self._ratfun.as_ZPK()
        result = K * scale
        for zero in zeros:
            result *= (1 - exp(zero * dt) / z)
        for pole in poles:
            result /= (1 - exp(pole * dt) / z)
        result *= undef

        result.is_causal = self.is_causal
        return result

    def impulse_invariance_transform(self):
        """This samples the impulse response and then calculates the
        Z-transform.  It does not work if the impulse response
        has Dirac deltas, say for a transfer function that is
        a pure delay or is not-strictly proper.

        The discrete-time and continuous-time impulse responses
        are identical at the sampling instants n * dt.

        The data needs to be sampled many times the bandwidth to avoid
        aliasing.

        See also `bilinear_transform` and `matched_ztransform`.

        See also `discretize` with regards to scaling of the result.
        """

        from .discretetime import n, z, dt

        # An alternative approach is to expand as partial fractions
        # and then replace (s + alpha) with (1 - exp(-alpha * dt) *
        # z**-1) in the denominator of each partial fraction (the
        # residues are unchanged).  This maps the pole at -alpha to
        # exp(-alpha * dt).

        h = self.ILT()
        if h.has(DiracDelta):
            raise ValueError('Impulse response has Dirac-deltas')

        scale = dt if self.is_ratio else 1

        hn = h.subs(n * dt)
        H = hn.ZT() * scale
        return H

    def discretize(self, method='bilinear', alpha=0.5):
        """Convert to a discrete-time approximation in the z-domain:

        :math:`H(z) \approx K H_c(s)`

        where :math:`K` is a scale factor.

        Note, the scaling is different for admittance, impedance, or
        transfer function expressions compared to voltage, current,
        and undefined expressions.

        The scaling is chosen so that the discrete-time voltage and
        current expressions have plots similar to the continuous-time
        forms.  For example, with the impulse-invariance method for a
        continuous-time voltage signal :math:`v_c(t)`:

        :math:`v[n] = v_c(n \Delta t)`

        However, for a transfer-function with a continuous-time
        impulse response :math:`h_c(t)`, then

        :math:`h[n] = \Delta t h_c(n \Delta t)`.

        This corrects the scaling when approximating a continuous-time
        convolution by a discrete-time convolution.

        The default method is 'bilinear'.  Other methods are:
        'impulse-invariance' 'bilinear', 'tustin', 'trapezoidal'
        'generalized-bilinear', 'gbf' controlled by the parameter
        `alpha` 'euler', 'forward-diff', 'forward-euler'
        'backward-diff', 'backward-euler' 'simpson', 'matched-Z',
        'zero-pole-matching'

        """

        if method in ('gbf', 'generalized-bilinear'):
            return self.generalized_bilinear_transform(alpha)
        elif method in ('bilinear', 'tustin', 'trapezoidal'):
            return self.generalized_bilinear_transform(0.5)
        elif method in ('euler', 'forward-diff', 'forward-euler'):
            return self.generalized_bilinear_transform(0)
        elif method in ('backward-diff', 'backward-euler'):
            return self.generalized_bilinear_transform(1)
        elif method in ('simpson', ):
            return self.simpson_transform()
        elif method in ('impulse-invariance', ):
            return self.impulse_invariance_transform()
        elif method in ('matched-Z', 'zero-pole-matching'):
            return self.matched_ztransform()
        else:
            raise ValueError('Unsupported method %s' % method)

    def zdomain(self, **assumptions):
        return self.discretize(**assumptions)

    def discrete_frequency(self, method='bilinear', **assumptions):
        return self.zdomain(method=method).discrete_frequency(**assumptions)

    def discrete_time(self, method='bilinear', **assumptions):
        return self.zdomain(method=method).discrete_time(**assumptions)


def tf(numer, denom=1, var=None):
    """Create a transfer function from lists of the coefficient
    for the numerator and denominator."""

    if var is None:
        var = ssym

    N = Poly(sympify(numer), var)
    D = Poly(sympify(denom), var)

    return LaplaceDomainTransferFunction(N / D, causal=True)


def zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists (or dictionaries) of zeros and poles,
    and from a constant gain."""

    if var is None:
        var = ssym
    return LaplaceDomainTransferFunction(_zp2tf(sympify(zeros), sympify(poles),
                                                sympify(K), var), causal=True)


def pr2tf(poles, residues, var=None):
    """Create a transfer function from lists of poles and residues."""

    if var is None:
        var = ssym
    return LaplaceDomainTransferFunction(_pr2tf(sympify(poles),
                                                sympify(residues), var),
                                         causal=True)


def sexpr(arg, **assumptions):
    """Create LaplaceDomainExpression object.  If `arg` is ssym return s"""

    if arg is ssym:
        return s
    return expr_make('laplace', arg, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

classes = expressionclasses.register('laplace', LaplaceDomainExpression)
LaplaceDomainVoltage = classes['voltage']
LaplaceDomainCurrent = classes['current']
LaplaceDomainAdmittance = classes['admittance']
LaplaceDomainImpedance = classes['impedance']
LaplaceDomainTransferFunction = classes['transfer']

from .texpr import TimeDomainExpression, texpr  # nopep8
from .phasor import PhasorDomainExpression  # nopep8
from .phasor import PhasorRatioDomainExpression  # nopep8

s = LaplaceDomainExpression('s')
s.units = uu.rad / uu.s

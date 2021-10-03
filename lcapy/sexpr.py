"""This module provides the LaplaceDomainExpression class to represent
s-domain (Laplace domain) expressions.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import LaplaceDomain
from .inverse_laplace import inverse_laplace_transform
from .sym import ssym, tsym, j, pi, sympify
from .ratfun import _zp2tf, _pr2tf, Ratfun
from .expr import Expr, symbol, expr, ExprDict, exprcontainer, expr_make
from .units import u as uu
from .functions import sqrt
import numpy as np
from sympy import limit, exp, Poly, Integral, div, oo, Eq, Expr as symExpr


__all__ = ('sexpr', 'zp2tf', 'tf', 'pr2tf')


class LaplaceDomainExpression(LaplaceDomain, Expr):
    """s-domain expression or symbol."""

    var = ssym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)                
        super(LaplaceDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr        
        if check and expr.has(tsym) and not expr.has(Integral):
            raise ValueError(
                's-domain expression %s cannot depend on t' % expr)

    def as_expr(self):
        return LaplaceDomainExpression(self)

    @classmethod
    def from_poles_residues(cls, poles, residues):
        """Create a transfer function from lists of poles and residues.

        See also from_zeros_poles_gain, from_numer_denom"""        

        return cls(pr2tf(poles, residues, cls.var), causal=True)

    @classmethod
    def from_zeros_poles_gain(cls, zeros, poles, K=1):
        """Create a transfer function from lists of zeros and poles,
        and from a constant gain.

        See also from_poles_residues, from_numer_denom"""        

        return cls(zp2tf(zeros, poles, K, cls.var), causal=True)

    @classmethod
    def from_numer_denom(cls, numer, denom):
        """Create a transfer function from lists of the coefficient
        for the numerator and denominator.

        See also from_zeros_poles_gain, from_poles_residues"""        

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

    def inverse_laplace(self, **assumptions):
        """Attempt inverse Laplace transform.

        If causal=True the response is zero for t < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for t < 0.
        Otherwise the result is only known for t >= 0.

        """
        assumptions = self.assumptions.merge(**assumptions)
        result = inverse_laplace_transform(self.expr, self.var, tsym,
                                           **assumptions)

        return self.change(result, domain='time', units_scale=uu.Hz, **assumptions)

    def ILT(self, **assumptions):
        """Attempt inverse Laplace transform.

        If causal=True the response is zero for t < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for t < 0.
        Otherwise the result is only known for t >= 0.

        """

        return self.inverse_laplace(**assumptions)
        
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

        if self.is_causal or assumptions.get('causal', False):
            # Note, this does not apply for 1 / s.            
            tmp = self(jw)
            if tmp.real != 0:
                return self.change(tmp(2 * pi * f), domain='fourier',
                                   **assumptions)                        
        
        result = self.time(**assumptions).fourier(**assumptions)
        return result

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import jw
        
        if self.is_causal:
            # Note, this does not apply for 1 / s.
            tmp = self(jw)
            if tmp.real != 0:
                return self.change(tmp, domain='angular fourier',
                                   **assumptions)                
        
        result = self.time(**assumptions).angular_fourier(**assumptions)
        return result

    def norm_angular_fourier(self, **assumptions):
        """Convert to normalized angular Fourier domain."""
        from .symbols import jw, Omega
        from .dsym import dt
        
        if self.is_causal:
            # Note, this does not apply for 1 / s.
            tmp = self(j * Omega / dt)
            if tmp.real != 0:
                return self.change(tmp, domain='norm angular fourier',
                                   **assumptions)                
        
        result = self.time(**assumptions).norm_angular_fourier(**assumptions)
        return result

    def norm_fourier(self, **assumptions):
        """Convert to normalized Fourier domain."""
        from .symbols import jw, F
        from .dsym import dt
        
        if self.is_causal:
            # Note, this does not apply for 1 / s.
            tmp = self(j * F / dt)
            if tmp.real != 0:
                return self.change(tmp, domain='norm fourier',
                                   **assumptions)                
        
        result = self.time(**assumptions).norm_fourier(**assumptions)
        return result        

    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        result = PhasorFrequencyDomainExpression.from_laplace(self, **assumptions)
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

    def angular_frequency_response(self, wvector=None):
        """Convert to angular frequency domain and evaluate response if
        angular frequency vector specified.

        """
        from .symbols import omega        

        X = self.subs(j * omega)

        if wvector is None:
            return X

        return X.evaluate(wvector)

    def frequency_response(self, fvector=None):
        """Convert to frequency domain and evaluate response if frequency
        vector specified.

        """
        from .symbols import f        

        X = self.subs(j * 2 * pi * f)

        if fvector is None:
            return X

        return X.evaluate(fvector)

    def _response_adhoc(self, xvector, tvector):    
        """Evaluate response of system with applied signal.
        This method assumes that the data is well over-sampled."""
        
        # Perform polynomial long division so expr1 = Q + M / D                
        N, D, delay = self._decompose()
        Q, M = div(N, D)
        expr1 = M / D

        Nt = len(tvector)

        # Evaluate impulse response.
        ndt = tvector[1] - tvector[0]        
        th = np.arange(Nt) * ndt
        H= LaplaceDomainExpression(expr1, **self.assumptions)
        hvector = H.transient_response(th)

        ty = tvector
        y = np.convolve(xvector, hvector)[0:Nt] * ndt

        if Q:
            # Handle Dirac deltas and their derivatives.
            C = expr(Q).coeffs()
            for n, c in enumerate(reversed(C)):

                y += float(c.expr) * xvector

                xvector = np.diff(xvector) / ndt
                xvector = np.hstack((xvector, 0))

        from scipy.interpolate import interp1d

        if delay != 0.0:
            # Try linear interpolation; should oversample first...
            y = interp1d(ty, y, bounds_error=False, fill_value=0)
            y = y(tvector - delay)

        return y

    def _response_bilinear(self, xvector, tvector):

        from .dsym import dt
        import scipy.signal as signal        
        
        # TODO: fix time offset
        if tvector[0] != 0:
            raise ValueError('Need initial time of 0')
        
        ndt = tvector[1] - tvector[0]
        
        Hz = self.bilinear_transform().subs(dt, ndt)
        fil = Hz.dlti_filter()

        a = [a1.fval for a1 in fil.a]        
        b = [b1.fval for b1 in fil.b]
        
        y = signal.lfilter(b, a, xvector)
        return y

    def response(self, xvector, tvector, method='bilinear'):
        """Evaluate response to input signal `xvector` at times 
        `tvector`.  This returns a NumPy array."""

        symbols = self.symbols
        symbols.pop('s', None)
        if symbols != {}:
            raise ValueError('Have undefined symbols: %s' % symbols)

        if len(xvector) != len(tvector):
            raise ValueError('x must have same length as t')

        dt = tvector[1] - tvector[0]
        if not np.allclose(np.diff(tvector), np.ones(len(tvector) - 1) * dt):
            raise (ValueError, 't values not equally spaced')

        if self.is_constant:
            return float(self.expr) * xvector

        if method == 'adhoc':
            return self._response_adhoc(xvector, tvector)
        elif method == 'bilinear':
            return self._response_bilinear(xvector, tvector)
        raise ValueError('Unknown method %s' % method)
        
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

    def _decompose(self):

        N, D, delay = self._ratfun.as_ratfun_delay()                

        return N, D, delay

    def differential_equation(self, input='x', output='y'):
        """Create differential equation from transfer function. 

        For example,  
        >>> H = (s + 3) / (s**2 + 4)  
        >>> H.differential_equation()
                 d                    d       
        3.y(t) + --(y(t)) = 4.x(t) + ---(x(t))
                 dt                    2      
                                     dt       
        """

        H = self
        x = texpr('%s(t)' % input)
        y = texpr('%s(t)' % output)

        X = x.LT()
        Y = y.LT()

        N = self.N
        D = self.D
        
        lhs = (N * Y).ILT(causal=True)
        rhs = (D * X).ILT(causal=True)

        return TimeDomainExpression(Eq(lhs.expr, rhs.expr))

    def dlti_filter(self, method='bilinear'):
        """Create DLTI filter using bilinear transform."""

        if method != 'bilinear':
            raise ValueError('Unsupported transform ' + method)
        return self.bilinear_transform().simplify().dlti_filter()
    
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

    def bode_plot(self, fvector=None, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Bode
        plot (but without the straight line approximations).  fvector
        specifies the frequencies.  If it is a tuple (f1, f2), it sets
        the frequency limits.   Since a logarithmic frequency scale is used,
        f1 must be greater than 0.

        This method makes the assumption that the expression is causal.

        """        

        return self.fourier(causal=True).bode_plot(fvector, **kwargs)

    def nyquist_plot(self, fvector=None, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Nyquist
        plot (imaginary part versus real part).  fvector specifies the
        frequencies.  If it is a tuple (f1, f2), it sets the frequency
        limits.

        `npoints` set the number of plotted points.

        The unit circle is shown by default.  This can be disabled with `unitcircle=False`.

        This method makes the assumption that the expression is causal.

        """        

        return self.fourier(causal=True).nyquist_plot(fvector, **kwargs)

    def nichols_plot(self, fvector=None, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Nichols
        plot (dB versus phase).  fvector specifies the frequencies.
        If it is a tuple (f1, f2), it sets the frequency limits.

        `npoints` set the number of plotted points.

        This method makes the assumption that the expression is causal.

        """        

        return self.fourier(causal=True).nichols_plot(fvector, **kwargs)        

    def generalized_bilinear_transform(self, alpha=0.5):
        
        from .discretetime import z, dt

        if alpha < 0 or alpha > 1:
            raise ValueError("alpha must be between 0 and 1 inclusive")
        
        return self.subs((1 / dt) * (1 - z**-1) / (alpha + (1 - alpha) * z**-1))
        
    def bilinear_transform(self):
        """Approximate s = ln(z) / dt

        by s = (2 / dt) * (1 - z**-1) / (1 + z**-1)

        This is also called Tustin's method and is equivalent to the
        trapezoidal method."""

        # TODO: add frequency warping as an option
        return self.generalized_bilinear_transform(0.5)

    def forward_euler_transform(self):
        """Approximate s = ln(z)

        by s = (1 / dt) * (1 - z**-1) / z**-1"""

        return self.generalized_bilinear_transform(0)

    def backward_euler_transform(self):
        """Approximate s = ln(z)

        by s = (1 / dt) * (1 - z**-1)"""

        return self.generalized_bilinear_transform(1)

    def discretize(self, method='bilinear', alpha=0.5):
        """Convert to a discrete-time approximation.

        The default method is 'bilinear'.  Other methods are
        'forward_euler', 'backward_euler', and 'gbf'.
        The latter has a parameter `alpha`."""

        if method == 'gbf':
            return self.generalized_bilinear_transform(alpha)
        elif method in ('bilinear', 'tustin'):
            return self.generalized_bilinear_transform(0.5)
        elif method in ('euler', 'forward_diff', 'forward_euler'):
            return self.generalized_bilinear_transform(0)
        elif method in ('backward_diff', 'backward_euler'):
            return self.generalized_bilinear_transform(1)
        else:
            raise ValueError('Unsupported method %s' % method)        

 
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
    return LaplaceDomainTransferFunction(_pr2tf(sympify(poles), sympify(residues), var),
                                         causal=True)


def sexpr(arg, **assumptions):
    """Create LaplaceDomainExpression object.  If `arg` is ssym return s"""

    if arg is ssym:
        return s
    return expr_make('laplace', arg, **assumptions)


from .expressionclasses import expressionclasses

classes = expressionclasses.register('laplace', LaplaceDomainExpression)
LaplaceDomainVoltage = classes['voltage']
LaplaceDomainCurrent = classes['current']
LaplaceDomainAdmittance = classes['admittance']
LaplaceDomainImpedance = classes['impedance']
LaplaceDomainTransferFunction = classes['transfer']

from .texpr import TimeDomainExpression, texpr
from .phasor import PhasorFrequencyDomainExpression

s = LaplaceDomainExpression('s')
s.units = uu.rad / uu.s

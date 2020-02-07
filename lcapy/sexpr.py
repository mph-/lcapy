"""This file provides the sExpr class to represent s-domain (Laplace
domain) expressions.

Copyright 2014--2019 Michael Hayes, UCECE

"""

from __future__ import division
from .laplace import inverse_laplace_transform
from .sym import ssym, tsym, j, pi
from .vector import Vector
from .ratfun import _zp2tf, Ratfun
from .expr import Expr, symbol, expr, ExprDict
from .functions import sqrt
import sympy as sym
import numpy as np


__all__ = ('Hs', 'Is', 'Vs', 'Ys', 'Zs', 'zp2tf', 'tf')


class sExpr(Expr):
    """s-domain expression or symbol."""

    var = ssym

    def __init__(self, val, **assumptions):

        super(sExpr, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = tExpr

        if self.expr.find(tsym) != set():
            raise ValueError(
                's-domain expression %s cannot depend on t' % self.expr)

    def tdifferentiate(self):
        """Differentiate in t-domain (multiply by s)."""

        return self.__class__(self.expr * self.var)

    def tintegrate(self):
        """Integrate in t-domain (divide by s)."""

        return self.__class__(self.expr / self.var)

    def delay(self, T):
        """Apply delay of T seconds by multiplying by exp(-s T)."""

        T = self.__class__(T)
        return self.__class__(self.expr * sym.exp(-s * T))

    @property
    def jomega(self):
        """Return expression with s = j omega."""

        from .symbols import jomega
        return self.subs(self.var, jomega)

    def initial_value(self):
        """Determine value at t = 0."""

        return self.__class__(sym.limit(self.expr * self.var, self.var, sym.oo))

    def final_value(self):
        """Determine value at t = oo."""

        return self.__class__(sym.limit(self.expr * self.var, self.var, 0))

    def inverse_laplace(self, **assumptions):
        """Attempt inverse Laplace transform.

        If causal=True the response is zero for t < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for t < 0.
        Otherwise the result is only known for t >= 0.

        """

        if assumptions == {}:
            assumptions = self.assumptions.copy()

        result = inverse_laplace_transform(self.expr, self.var, tsym, **assumptions)

        if hasattr(self, '_laplace_conjugate_class'):
            result = self._laplace_conjugate_class(result)
        else:
            result = tExpr(result)
        return result

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.inverse_laplace(**assumptions)

    def laplace(self, **assumptions):
        """Convert to s-domain."""

        if assumptions == {}:
            assumptions = self.assumptions.copy()
        
        return self.__class__(self, **assumptions)

    def fourier(self, **assumptions):
        """Convert to Fourier domain."""
        from .symbols import f
        
        if assumptions.get('causal', self.is_causal):
            return self.subs(j * 2 * pi * f)

        return self.time(**assumptions).fourier(**assumptions)

    def phasor(self, **assumptions):

        return self.time(**assumptions).phasor(**assumptions)

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

    def response(self, x, t):
        """Evaluate response to input signal x at times t."""

        if len(x) != len(t):
            raise ValueError('x must have same length as t')

        dt = t[1] - t[0]
        if not np.allclose(np.diff(t), np.ones(len(t) - 1) * dt):
            raise (ValueError, 't values not equally spaced')

        # Perform polynomial long division so expr = Q + M / D                
        N, D, delay = self._decompose()
        Q, M = sym.div(N, D)
        expr = M / D

        N = len(t)

        # Evaluate transient response.
        th = np.arange(N) * dt - dt
        h = sExpr(expr).transient_response(th)

        print('Convolving...')
        ty = t
        y = np.convolve(x, h)[0:N] * dt

        if Q:
            # Handle Dirac deltas and their derivatives.
            C = Q.all_coeffs()
            for n, c in enumerate(C):

                y += c * x

                x = np.diff(x) / dt
                x = np.hstack((x, 0))

        from scipy.interpolate import interp1d

        if delay != 0.0:
            print('Interpolating...')
            # Try linear interpolation; should oversample first...
            y = interp1d(ty, y, bounds_error=False, fill_value=0)
            y = y(t - delay)

        return y

    def _decompose(self):

        N, D, delay = Ratfun(self, s).as_ratfun_delay()                

        return N, D, delay

    def evaluate(self, svector=None):

        return super(sExpr, self).evaluate(svector)

    def plot(self, t=None, **kwargs):
        """Plot pole-zero map."""

        from .plot import plot_pole_zero
        return plot_pole_zero(self, **kwargs)

    def parameterize(self, zeta=True):
        """Parameterize first and second-order expressions.

        For first order systems, parameterize as:

        K * (s + beta) / (s + alpha)

        K / (s + alpha)

        K (s + beta)

        where appropriate.

        If `zeta` is True, parameterize second-order expression in
        standard form using damping factor and natural frequency
        representation, i.e.

        N(s) / (s**2 + 2 * zeta * omega_0 * s + omega_0**2)
        
        otherwise parameterize as
        
        N(s) / (s**2 + 2 * sigma_1 * s + omega_1**2 + sigma_1**2)

        """

        def def1(defs, symbolname, value):
            sym1 = symbol(symbolname)
            defs[symbolname] = value
            return sym1

        ndegree = self.N.degree        
        ddegree = self.D.degree
        ncoeffs = self.N.coeffs(norm=True)
        dcoeffs = self.D.coeffs(norm=True)

        defs = ExprDict()

        K = self.K
        if ndegree < 1 and ddegree < 1:
            return self, defs
        if ndegree == 1 and ddegree == 1:
            K = def1(defs, 'K', K)
            alpha = def1(defs, 'alpha', dcoeffs[1])
            beta = def1(defs, 'beta', ncoeffs[1])
            return K * (s + beta) / (s + alpha), defs
        if ndegree == 1 and ddegree == 0:
            K = def1(defs, 'K', K)
            beta = def1(defs, 'beta', ncoeffs[1])
            return K * (s + beta), defs
        if ndegree == 0 and ddegree == 1:
            K = def1(defs, 'K', K)
            alpha = def1(defs, 'alpha', dcoeffs[1])
            return K / (s + alpha), defs
        if ddegree == 2:
            K = def1(defs, 'K', K)
            coeffs = self.N.coeffs()

            if not zeta:
                sigma1 = def1(defs, 'sigma_1', dcoeffs[1] / 2)
                omega1 = def1(defs, 'omega_1',
                              sqrt(dcoeffs[2] - (dcoeffs[1] / 2)**2).simplify())
                return K * (self.N / coeffs[0]) / (s**2 + 2 * sigma1 * s + sigma1**2 + omega1**2), defs
                
            omega0 = def1(defs, 'omega_0', sqrt(dcoeffs[2]))
            zeta = def1(defs, 'zeta', dcoeffs[1] / (2 * sqrt(dcoeffs[2])))
            return K * (self.N / coeffs[0]) / (s**2 + 2 * zeta * omega0 * s + omega0**2), defs
        
        return self, defs

    def bilinear_transform(self):

        from .discretetime import z, dt

        # s = ln(z) / dt gives the exact solution
        
        return self.subs((2 / dt) * (1 - z**-1) / (1 + z**-1))
        
    
    
# Perhaps use a factory to create the following classes?

class Zs(sExpr):

    """s-domain impedance value."""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, causal=True, **assumptions):

        super(Zs, self).__init__(val, causal=causal, **assumptions)
        self._laplace_conjugate_class = Zt

    def cpt(self):
        from .oneport import R, C, L, Z

        if self.is_number or self.is_dc:
            return R(self.expr)

        z = self * s

        if z.is_number:
            return C((1 / z).expr)

        z = self / s

        if z.is_number:
            return L(z.expr)

        return Z(self)


class Ys(sExpr):

    """s-domain admittance value."""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, causal=True, **assumptions):

        super(Ys, self).__init__(val, causal=causal, **assumptions)
        self._laplace_conjugate_class = Yt

    def cpt(self):
        from .oneport import G, C, L, Y

        if self.is_number or self.is_dc:
            return G(self.expr)

        y = self * s

        if y.is_number:
            return L((1 / y).expr)

        y = self / s

        if y.is_number:
            return C(y.expr)

        return Y(self)


class Vs(sExpr):

    """s-domain voltage (units V s / radian)."""

    quantity = 's-Voltage'
    units = 'V/Hz'

    def __init__(self, val, **assumptions):

        super(Vs, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Vt

    def cpt(self):
        from .oneport import V
        return V(self)


class Is(sExpr):

    """s-domain current (units A s / radian)."""

    quantity = 's-Current'
    units = 'A/Hz'

    def __init__(self, val, **assumptions):

        super(Is, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = It

    def cpt(self):
        from .oneport import I
        
        return I(self)


class Hs(sExpr):

    """s-domain ratio"""

    quantity = 's-ratio'
    units = ''

    def __init__(self, val, **assumptions):

        super(Hs, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Ht

        
class VsVector(Vector):

    _typewrap = Vs


class IsVector(Vector):

    _typewrap = Is


class YsVector(Vector):

    _typewrap = Ys


class ZsVector(Vector):

    _typewrap = Zs


def tf(numer, denom=1, var=None):
    """Create a transfer function from lists of the coefficient
    for the numerator and denominator."""

    if var is None:
        var = ssym

    N = sym.Poly(numer, var)
    D = sym.Poly(denom, var)

    return Hs(N / D)


def zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists of zeros and poles,
    and from a constant gain."""

    if var is None:
        var = ssym
    return Hs(_zp2tf(zeros, poles, K, var))


def sexpr(arg):
    """Create sExpr object.  If `arg` is ssym return s"""

    if arg is ssym:
        return s
    return sExpr(arg)


from .texpr import Ht, It, Vt, Yt, Zt, tExpr
s = sExpr('s')


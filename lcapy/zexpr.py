"""This file provides the zExpr class to represent z-domain expressions.

Copyright 2020 Michael Hayes, UCECE

"""

from __future__ import division
#from .ztransform import inverse_ztransform
from .sym import j, pi
from .dsym import nsym, ksym, zsym, dt
from .vector import Vector
from .ratfun import _zp2tf, Ratfun
from .expr import Expr, symbol, expr, ExprDict
from .functions import sqrt, exp
import sympy as sym
import numpy as np

__all__ = ('Hz', 'Iz', 'Vz', 'Yz', 'Zz')

class zExpr(Expr):
    """z-domain expression or symbol."""

    var = zsym

    def __init__(self, val, **assumptions):

        super(zExpr, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = None

        if self.expr.find(nsym) != set():
            raise ValueError(
                'z-domain expression %s cannot depend on n' % self.expr)

    def differentiate(self):
        """Differentiate (multiply by z)."""

        return self.__class__(self.expr * self.var)

    def integrate(self):
        """Integrate (divide by z)."""

        return self.__class__(self.expr / self.var)

    def delay(self, T):
        """Apply delay of T seconds by multiplying by exp(-s T)."""

        T = self.__class__(T)
        return self.__class__(self.expr * sym.exp(-s * T))

    def initial_value(self):
        """Determine value at n = 0."""

        return self.__class__(sym.limit(self.expr * self.var, self.var, sym.oo))

    def final_value(self):
        """Determine value at n = oo."""

        return self.__class__(sym.limit(self.expr * self.var, self.var, 0))

    def inverse_ztransform(self, **assumptions):
        """Attempt inverse Znransform transform.

        If causal=True the response is zero for n < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for n < 0.
        Otherwise the result is only known for n >= 0.

        """

        if assumptions == {}:
            assumptions = self.assumptions.copy()

        result = inverse_ztransform_transform(self.expr, self.var, tsym, **assumptions)

        if hasattr(self, '_ztransform_conjugate_class'):
            result = self._ztransform_conjugate_class(result)
        else:
            result = nexpr(result)
        return result

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return self.inverse_ztransform(**assumptions)

    def ztransform(self, **assumptions):
        """Convert to z-domain."""

        if assumptions == {}:
            assumptions = self.assumptions.copy()
        
        return self.__class__(self, **assumptions)

    def DTFT(self, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""
        from .symbols import f
        
        if assumptions.get('causal', self.is_causal):
            return self.subs(exp(j * 2 * pi * f * dt))

        raise RuntimeError('TODO')
        #return self.discrete_time(**assumptions).discrete_fourier(**assumptions)

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
        h = zExpr(expr).transient_response(th)

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

        return super(zExpr, self).evaluate(svector)

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
    
    def inverse_bilinear_transform(self):

        from .symbols import s
        from .discretetime import dt

        # z = exp(s * dt) gives the exact solution
        
        return self.subs((1 + s * dt / 2) / (1 - s * dt / 2))

    
# Perhaps use a factory to create the following classes?

class Zz(zExpr):

    """z-domain impedance value."""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, causal=True, **assumptions):

        super(Zz, self).__init__(val, causal=causal, **assumptions)
        self._ztransform_conjugate_class = Zn


class Yz(zExpr):

    """z-domain admittance value."""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, causal=True, **assumptions):

        super(Yz, self).__init__(val, causal=causal, **assumptions)
        self._ztransform_conjugate_class = Yn


class Vz(zExpr):

    """z-domain voltage (units V s / radian)."""

    quantity = 'z-Voltage'
    units = 'V/Hz'

    def __init__(self, val, **assumptions):

        super(Vz, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Vn


class Iz(zExpr):

    """z-domain current (units A s / radian)."""

    quantity = 'z-Current'
    units = 'A/Hz'

    def __init__(self, val, **assumptions):

        super(Iz, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = In


class Hz(zExpr):

    """z-domain ratio"""

    quantity = 'z-ratio'
    units = ''

    def __init__(self, val, **assumptions):

        super(Hz, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Hn

        
class VzVector(Vector):

    _typewrap = Vz


class IzVector(Vector):

    _typewrap = Iz


class YzVector(Vector):

    _typewrap = Yz


class ZzVector(Vector):

    _typewrap = Zz


def zexpr(arg):
    """Create zExpr object.  If `arg` is zsym return z"""

    if arg is zsym:
        return z
    return zExpr(arg)


from .nexpr import Hn, In, Vn, Yn, Zn, nExpr
z = zExpr('z')


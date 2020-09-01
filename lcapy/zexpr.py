"""This module provides the zExpr class to represent z-domain expressions.

Copyright 2020 Michael Hayes, UCECE

"""

from __future__ import division
from .ztransform import inverse_ztransform
from .sym import j, pi
from .dsym import nsym, ksym, zsym, dt
from .vector import Vector
from .ratfun import _zp2tf, Ratfun
from .dexpr import dExpr
from .expr import symbol, expr, ExprDict
from .functions import sqrt, exp
import numpy as np
from sympy import Eq, div, limit, oo, Sum


__all__ = ('Hz', 'Iz', 'Vz', 'Yz', 'Zz')

class zExpr(dExpr):
    """z-domain expression or symbol."""

    var = zsym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)

        super(zExpr, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = nExpr

        expr = self.expr
        if check and expr.find(nsym) != set() and not expr.has(Sum):
            raise ValueError(
                'z-domain expression %s cannot depend on n' % expr)
        if check and expr.find(ksym) != set() and not expr.has(Sum):
            raise ValueError(
                'z-domain expression %s cannot depend on k' % expr)

    def ndifferentiate(self):
        """First order difference in n-domain."""

        q = 1 / (1 - 1 / self.var) * dt
        return self.__class__(self.expr / q, **self.assumptions)

    def nintegrate(self):
        """First order integration in n-domain."""

        q = 1 / (1 - 1 / self.var) * dt       
        return self.__class__(self.expr * q, **self.assumptions)

    def initial_value(self):
        """Determine value at n = 0."""

        return self.__class__(limit(self.expr * self.var, self.var, oo))

    def final_value(self):
        """Determine value at n = oo."""

        return self.__class__(limit(self.expr * self.var, self.var, 0))

    def inverse_ztransform(self, **assumptions):
        """Attempt inverse Z ransform.

        If causal=True the response is zero for n < 0 and
        the result is multiplied by UnitStep(n)
        If ac=True or dc=True the result is extrapolated for n < 0.
        Otherwise the result is only known for n >= 0.

        """

        if assumptions == {}:
            assumptions = self.assumptions.copy()

        result = inverse_ztransform(self.expr, self.var, nsym, **assumptions)

        if hasattr(self, '_ztransform_conjugate_class'):
            result = self._ztransform_conjugate_class(result)
        else:
            result = nexpr(result)
        return result

    def IZT(self, **assumptions):
        return self.inverse_ztransform(**assumptions)

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

        q = 1 / (1 - 1 / self.var)
        H = self.__class__(self * q, **self.assumptions)
        return H.transient_response(tvector)

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
        Q, M = div(N, D)
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

        if 'unitcircle' not in kwargs:
            kwargs['unitcircle'] = True

        from .plot import plot_pole_zero
        return plot_pole_zero(self, **kwargs)

    def inverse_bilinear_transform(self):

        from .symbols import s
        from .discretetime import dt

        # z = exp(s * dt) gives the exact solution
        
        return self.subs((1 + s * dt / 2) / (1 - s * dt / 2))

    def discrete_time_fourier_transform(self, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""

        from .symbols import f

        if assumptions.get('causal', self.is_causal):
            return self.subs(exp(j * 2 * pi * f * dt))

        return self.IZT(**assumptions).DTFT()

    def DTFT(self, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""
    
        return self.discrete_time_fourier_transform(**assumptions) 

    def decompose_AB(self):

        C, R = self.factor_const()
        
        invz = expr('invz')
        H = R.replace(z, 1 / invz).factor()
        r = Ratfun(H, invz.expr)
        B = zExpr(r.N).replace(invz, 1 / z)
        A = zExpr(r.D).replace(invz, 1 / z)

        C1, R1 = A.term_const()
        if C1.is_negative:
            A = -A
            B = -B
        
        return A, B * C
    
    def difference_equation(self, input='x', output='y', form='iir'):
        """Create difference equation from transfer function.

        form can be 'fir' or 'iir' ('direct form I').
        """

        H = self
        x = nexpr('%s(n)' % input)
        y = nexpr('%s(n)' % output)

        X = x.ZT()
        Y = y.ZT()

        if form in ('iir', 'direct form I'):
            # Direct form I
            A, B = self.decompose_AB()
            
            lhs = (A * Y).IZT(causal=True)
            rhs = (B * X).IZT(causal=True)
            
        elif form == 'fir':
            H = H.partfrac()
            lhs = y
            rhs = (H * X).IZT(causal=True)

        else:
            raise ValueError('Unhandled form ' + form)    

        return nExpr(Eq(lhs.expr, rhs.expr))        
        
    
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


def zexpr(arg, **assumptions):
    """Create zExpr object.  If `arg` is zsym return z"""

    if arg is zsym:
        return z
    return zExpr(arg, **assumptions)


from .nexpr import Hn, In, Vn, Yn, Zn, nExpr, nexpr
z = zExpr('z')

"""This module provides the ZDomainExpression class to represent z-domain expressions.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import ZDomain
from .ztransform import inverse_ztransform
from .sym import j, pi
from .dsym import nsym, ksym, zsym, dt
from .vector import Vector
from .ratfun import _zp2tf, Ratfun
from .dexpr import DiscreteExpression
from .expr import symbol, expr, ExprDict
from .functions import sqrt, exp
import numpy as np
from sympy import Eq, div, limit, oo, Sum


__all__ = ('zexpr', )

class ZDomainExpression(ZDomain, DiscreteExpression):
    """z-domain expression or symbol."""

    var = zsym


    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)

        super(ZDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr
        if check and expr.find(nsym) != set() and not expr.has(Sum):
            raise ValueError(
                'z-domain expression %s cannot depend on n' % expr)
        if check and expr.find(ksym) != set() and not expr.has(Sum):
            raise ValueError(
                'z-domain expression %s cannot depend on k' % expr)

    def as_expr(self):
        return ZDomainExpression(self)

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

        assumptions = self.assumptions.merge(**assumptions)
        result = inverse_ztransform(self.expr, self.var, nsym, **assumptions)
        return self.change(result, domain='discrete time', **assumptions)

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
        h = ZDomainExpression(expr).transient_response(th)

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

        N, D, delay = Ratfun(self, z).as_ratfun_delay()                

        return N, D, delay

    def evaluate(self, svector=None):

        return super(ZDomainExpression, self).evaluate(svector)

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
        B = ZDomainExpression(r.N).replace(invz, 1 / z)
        A = ZDomainExpression(r.D).replace(invz, 1 / z)

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

        return DiscreteTimeDomainExpression(Eq(lhs.expr, rhs.expr))        
        
    
def zexpr(arg, **assumptions):
    """Create ZDomainExpression object.  If `arg` is zsym return z"""

    if arg is zsym:
        return z
    return ZDomainExpression(arg, **assumptions)


from .expressionclasses import expressionclasses

expressionclasses.register('Z', ZDomainExpression)

from .nexpr import DiscreteTimeDomainExpression, nexpr
z = ZDomainExpression('z')

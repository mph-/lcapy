"""This module provides the ZDomainExpression class to represent z-domain expressions.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import ZDomain
from .inverse_ztransform import inverse_ztransform
from .sym import j, pi, fsym, omegasym
from .dsym import nsym, ksym, zsym, dt
from .vector import Vector
from .ratfun import _zp2tf, Ratfun
from .dexpr import DiscreteExpression
from .expr import symbol, expr, ExprDict
from .diffeq import DifferenceEquation
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
        """Plot pole-zero map.

        kwargs include:
        unitcircle - if True, draw unit circle
        axes - the plot axes to use otherwise a new figure is created
        xlabel - the x-axis label (default Re(z))
        ylabel - the y-axis label (default Im(z))
        xscale - the x-axis scaling
        yscale - the y-axis scaling
        in addition to those supported by the matplotlib plot command.
        
        The plot axes are returned."""

        if 'unitcircle' not in kwargs:
            kwargs['unitcircle'] = True

        from .plot import plot_pole_zero
        return plot_pole_zero(self, **kwargs)

    def inverse_bilinear_transform(self):

        from .symbols import s
        from .discretetime import dt

        # z = exp(s * dt) gives the exact solution
        
        return self.subs((1 + s * dt / 2) / (1 - s * dt / 2))

    def discrete_time_fourier_transform(self, norm=False, angular=False,
                                        **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""

        from .symbols import f

        if assumptions.get('causal', self.is_causal):
            result = self.subs(exp(j * 2 * pi * f * dt))
        else:
            return self.IZT(**assumptions).DTFT(norm=norm, angular=angular)

        if norm:
            result = result.subs(dt, 1)
        if angular:
            result = result.angular_fourier()            
        return result            

    def DTFT(self, norm=False, angular=False, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""
    
        return self.discrete_time_fourier_transform(**assumptions,
                                                    norm=norm, angular=angular)

    def as_ab(self):
        """Return lists of denominator and numerator coefficients
        when the denominator and numerator are expressed as polynomials
        in z**-1.  The lowest order coefficients are returned first."""

        C, R = self.factor_const()
        
        zi = symbol('zi')
        H = R.replace(z, 1 / zi).cancel()
        a = H.D.coeffs(zi)
        b = H.N.coeffs(zi)
        return a[::-1], list(np.array(b) * C)[::-1]

    def as_AB(self):

        C, R = self.factor_const()
        
        zi = symbol('zi')
        H = R.replace(z, 1 / zi).factor()
        r = Ratfun(H, zi.expr)
        B = ZDomainExpression(r.N).replace(zi, 1 / z)
        A = ZDomainExpression(r.D).replace(zi, 1 / z)

        C1, R1 = A.term_const()
        if C1.is_negative:
            A = -A
            B = -B
        
        return A, B * C
    
    def difference_equation(self, inputsym='x', outputsym='y', form='iir'):
        """Create difference equation from transfer function.

        `form` can be 'fir' or 'iir' ('direct form I').
        """

        H = self
        x = nexpr('%s(n)' % inputsym)
        y = nexpr('%s(n)' % outputsym)

        X = x.ZT()
        Y = y.ZT()

        if form in ('iir', 'direct form I'):
            # Direct form I
            return self.dlti_filter().difference_equation()
            
        elif form == 'fir':
            H = H.partfrac()
            lhs = y
            rhs = (H * X).IZT(causal=True)

        else:
            raise ValueError('Unhandled form ' + form)    

        return DifferenceEquation(lhs, rhs, inputsym, outputsym)

    def dlti_filter(self):
        """Create discrete-time linear time-invariant filter from discrete-time
        transfer function."""

        # TODO, perhaps add only to DiscreteTimeDomainTransfer?
        
        from .dltifilter import DLTIFilter
        
        if not self.is_rational_function:
            raise ValueError("Not a rational function")            
            
        N = self.N
        D = self.D
        n_n = N.coeffs()
        d_n = D.coeffs()
    
        if len(n_n) > len(d_n):
            raise ValueError("System not causal")

        bn = (len(d_n) - len(n_n)) * [0] + n_n
        an = d_n
        lpf = DLTIFilter(bn, an) 
        return lpf
        
def zexpr(arg, **assumptions):
    """Create ZDomainExpression object.  If `arg` is zsym return z"""

    from .expr import Expr
    
    if arg is zsym:
        return z

    if isinstance(arg, Expr):
        if assumptions == {}:
            return arg
        return arg.__class__(arg, **assumptions)
    
    return ZDomainExpression(arg, **assumptions)


from .expressionclasses import expressionclasses

expressionclasses.register('Z', ZDomainExpression)

from .nexpr import DiscreteTimeDomainExpression, nexpr
z = ZDomainExpression('z')

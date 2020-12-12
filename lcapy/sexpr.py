"""This module provides the LaplaceDomainExpression class to represent s-domain (Laplace
domain) expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .laplace import inverse_laplace_transform
from .sym import ssym, tsym, j, pi, sympify
from .ratfun import _zp2tf, _pr2tf, Ratfun
from .expr import Expr, symbol, expr, ExprDict, exprcontainer
from .functions import sqrt
import numpy as np
from sympy import limit, exp, Poly, Integral, div, oo, Eq, Expr as symExpr


__all__ = ('zp2tf', 'tf', 'pr2tf')


class LaplaceDomainExpression(Expr):
    """s-domain expression or symbol."""

    var = ssym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)                
        super(LaplaceDomainExpression, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = TimeDomainExpression

        expr = self.expr        
        if check and expr.find(tsym) != set() and not expr.has(Integral):
            raise ValueError(
                's-domain expression %s cannot depend on t' % expr)
        if check and expr.find(tsym) != set() and not expr.has(Integral):
            raise ValueError(
                's-domain expression %s cannot depend on t' % self.expr)

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

        assumptions = self.merge_assumptions(**assumptions)
        result = inverse_laplace_transform(self.expr, self.var, tsym,
                                           **assumptions)

        if hasattr(self, '_laplace_conjugate_class'):
            result = self._laplace_conjugate_class(result)
        else:
            result = TimeDomainExpression(result)
        return result

    def ILT(self, **assumptions):
        """Convert to t-domain.   This is an alias for inverse_laplace."""

        return self.inverse_laplace(**assumptions)
    
    def time(self, **assumptions):
        """Convert to time domain."""

        try:
            return self.inverse_laplace(**assumptions)
        except ValueError:
            return self.as_sum().inverse_laplace(**assumptions)            

    def laplace(self, **assumptions):
        """Convert to s-domain."""

        assumptions = self.merge_assumptions(**assumptions)
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
        Q, M = div(N, D)
        expr = M / D

        N = len(t)

        # Evaluate transient response.
        th = np.arange(N) * dt - dt
        h = LaplaceDomainExpression(expr).transient_response(th)

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

    def evaluate(self, svector=None):

        return super(LaplaceDomainExpression, self).evaluate(svector)

    def plot(self, t=None, **kwargs):
        """Plot pole-zero map."""

        from .plot import plot_pole_zero
        return plot_pole_zero(self, **kwargs)

    def parameterize(self, zeta=True):
        """Parameterize first and second-order expressions.

        For example, pexpr, defs = expr.parameterize()

        If parameterization is successful, defs is a dictionary
        of the paramters.  The original expression can be obtained
        with pexpr.subs(defs)

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
            from .cexpr import ConstantExpression
            
            sym1 = symbol(symbolname)
            defs[symbolname] = ConstantExpression(value)
            return sym1

        factors = self.as_ordered_factors()

        spowers = [s**-4, s**-3, s**-2, s**-1, s, s**2, s**3, s**4]
        for spower in spowers:
            if spower in factors:
                result, defs = (self / spower).parameterize(zeta)
                return result * spower, defs
        
        N = self.N
        D = self.D
        
        ndegree = N.degree        
        ddegree = D.degree
        ncoeffs = N.coeffs(norm=True)
        dcoeffs = D.coeffs(norm=True)

        result = None
        defs = ExprDict()

        K = self.K
        if ndegree < 1 and ddegree < 1:
            result = self
        elif ndegree == 1 and ddegree == 1:
            K = def1(defs, 'K', K)
            alpha = def1(defs, 'alpha', dcoeffs[1])
            beta = def1(defs, 'beta', ncoeffs[1])
            result = K * (s + beta) / (s + alpha)
        elif ndegree == 1 and ddegree == 0:
            K = def1(defs, 'K', K)
            beta = def1(defs, 'beta', ncoeffs[1])
            result = K * (s + beta)
        elif ndegree == 0 and ddegree == 1:
            K = def1(defs, 'K', K)
            alpha = def1(defs, 'alpha', dcoeffs[1])
            result = K / (s + alpha)
        elif ddegree == 2:
            K = def1(defs, 'K', K)
            coeffs = self.N.coeffs()

            if not zeta:
                sigma1 = def1(defs, 'sigma_1', dcoeffs[1] / 2)
                omega1 = def1(defs, 'omega_1',
                              sqrt(dcoeffs[2] - (dcoeffs[1] / 2)**2).simplify())
                result = K * (self.N / coeffs[0]) / (s**2 + 2 * sigma1 * s + sigma1**2 + omega1**2)
            else:
                omega0 = def1(defs, 'omega_0', sqrt(dcoeffs[2]))
                zeta = def1(defs, 'zeta', dcoeffs[1] / (2 * sqrt(dcoeffs[2])))
                result = K * (self.N / coeffs[0]) / (s**2 + 2 * zeta * omega0 * s + omega0**2)

        if result is None:
            # Copy?
            result = self

        return self.__class__(result, **self.assumptions), defs

    def bilinear_transform(self):
        """Approximate s = ln(z)

        by s = (2 / dt) * (1 - z**-1) / (1 + z**-1)

        This is also called Tustin's method and is equivalent to the
        trapezoidal method."""

        # TODO: add frequency warping as an option

        from .discretetime import z, dt

        return self.subs((2 / dt) * (1 - z**-1) / (1 + z**-1))

    def forward_euler_transform(self):
        """Approximate s = ln(z)

        by s = (1 / dt) * (1 - z**-1) / z**-1"""

        from .discretetime import z, dt
        
        return self.subs((1 / dt) * (1 - z**-1) / (z**-1))

    def backward_euler_transform(self):
        """Approximate s = ln(z)

        by s = (1 / dt) * (1 - z**-1)"""

        from .discretetime import z, dt
        
        return self.subs((1 / dt) * (1 - z**-1))        
    
    def transform(self, arg, **assumptions):
        """Transform into a different domain."""

        from .fexpr import FourierDomainExpression, f
        from .omegaexpr import AngularFourierDomainExpression, omega
        from .symbols import jomega        

        arg = expr(arg)

        is_causal = assumptions.get('causal', self.is_causal)
        
        if arg.has(jomega):
            if not is_causal:
                raise ValueError('Cannot convert non-causal s-expression to jomega domain')
            result = self.subs(jomega)
            # Handle args like 5 * jomega,  This might be too cute.
            return result.subs(arg / j, **assumptions)

        elif isinstance(arg, AngularFourierDomainExpression):
            if not is_causal:
                raise ValueError('Cannot convert non-causal s-expression to omega domain')
            result = self.subs(jomega)
            # Handle args like 5 * omega,  This might be too cute.
            return result.subs(arg, **assumptions)

        elif isinstance(arg, FourierDomainExpression):
            if not is_causal:
                raise ValueError('Cannot convert non-causal s-expression to f domain')
            result = self.subs(j * 2 * pi * f)
            # Handle args like 5 * f,  This might be too cute.
            return result.subs(arg, **assumptions)        
        
        return super(LaplaceDomainExpression, self).transform(arg, **assumptions)    

    
# Perhaps use a factory to create the following classes?

class LaplaceDomainImpedance(LaplaceDomainExpression):

    """s-domain impedance value."""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, causal=True, **assumptions):

        super(LaplaceDomainImpedance, self).__init__(val, causal=causal, **assumptions)
        self._laplace_conjugate_class = TimeDomainImpedance

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

    def network(self, form='default'):
        """Synthesise a network with an equivalent impedance.
        `form` includes: cauerI, cauerII, fosterI, fosterII.

        Note some methods generate networks with negative value
        components."""

        from .synthesis import network
        
        return network(self, form)


class LaplaceDomainAdmittance(LaplaceDomainExpression):

    """s-domain admittance value."""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, causal=True, **assumptions):

        super(LaplaceDomainAdmittance, self).__init__(val, causal=causal, **assumptions)
        self._laplace_conjugate_class = TimeDomainAdmittance

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

    def network(self, form='default'):
        """Synthesise a network with an equivalent impedance.
        `form` includes: cauerI, cauerII, fosterI, fosterII.

        Note some methods generate networks with negative value
        components."""        

        from .synthesis import network
        
        return network(1 / self, form)

    
class LaplaceDomainVoltage(LaplaceDomainExpression):

    """s-domain voltage (units V s / radian)."""

    quantity = 's-Voltage'
    units = 'V/Hz'
    superkind = 'Voltage'

    def __init__(self, val, **assumptions):

        super(LaplaceDomainVoltage, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = TimeDomainVoltage

    def cpt(self):
        from .oneport import V
        return V(self)

    def __mul__(self, x):
        """Multiply"""

        if isinstance(x, LaplaceDomainAdmittance):
            return LaplaceDomainCurrent(super(LaplaceDomainVoltage, self).__mul__(x))
        if isinstance(x, (ConstantExpression, LaplaceDomainExpression, symExpr, int, float, complex)):
            return super(LaplaceDomainVoltage, self).__mul__(x)
        self._incompatible(x, '*')

    def __truediv__(self, x):
        """Divide"""

        if isinstance(x, LaplaceDomainImpedance):
            return LaplaceDomainCurrent(super(LaplaceDomainVoltage, self).__truediv__(x))
        if isinstance(x, LaplaceDomainVoltage):
            return LaplaceDomainTransferFunction(super(LaplaceDomainVoltage, self).__truediv__(x))        
        if isinstance(x, (ConstantExpression, LaplaceDomainExpression,
                          symExpr, int, float, complex)):
            return super(LaplaceDomainVoltage, self).__truediv__(x)
        self._incompatible(x, '/')        

        
class LaplaceDomainCurrent(LaplaceDomainExpression):

    """s-domain current (units A s / radian)."""

    quantity = 's-Current'
    units = 'A/Hz'
    superkind = 'Current'

    def __init__(self, val, **assumptions):

        super(LaplaceDomainCurrent, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = TimeDomainCurrent

    def cpt(self):
        from .oneport import I
        
        return I(self)

    def __mul__(self, x):
        """Multiply"""

        if isinstance(x, LaplaceDomainImpedance):
            return LaplaceDomainVoltage(super(LaplaceDomainCurrent, self).__mul__(x))            
        if isinstance(x, (ConstantExpression, LaplaceDomainExpression, symExpr, int, float, complex)):
            return super(LaplaceDomainCurrent, self).__mul__(x)
        self._incompatible(x, '*')        

    def __truediv__(self, x):
        """Divide"""

        if isinstance(x, LaplaceDomainAdmittance):
            return LaplaceDomainVoltage(super(LaplaceDomainCurrent, self).__truediv__(x))
        if isinstance(x, LaplaceDomainCurrent):
            return LaplaceDomainTransferFunction(super(LaplaceDomainCurrent, self).__truediv__(x))                
        if isinstance(x, (ConstantExpression, LaplaceDomainExpression, symExpr, int, float, complex)):
            return super(LaplaceDomainCurrent, self).__truediv__(x)
        self._incompatible(x, '/')                
    

class LaplaceDomainTransferFunction(LaplaceDomainExpression):

    """s-domain ratio"""

    quantity = 's-ratio'
    units = ''

    def __init__(self, val, **assumptions):

        super(LaplaceDomainTransferFunction, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = TimeDomainImpulseResponse


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
    return LaplaceDomainTransferFunction(_pr2tf(sympify(poles), sympify(residues), var), causal=True)


def sexpr(arg, **assumptions):
    """Create LaplaceDomainExpression object.  If `arg` is ssym return s"""

    if arg is ssym:
        return s
    return LaplaceDomainExpression(arg, **assumptions)


from .texpr import TimeDomainImpulseResponse, TimeDomainCurrent
from .texpr import TimeDomainVoltage, TimeDomainAdmittance, TimeDomainImpedance
from .texpr import TimeDomainExpression, texpr
from .cexpr import ConstantExpression
s = LaplaceDomainExpression('s')


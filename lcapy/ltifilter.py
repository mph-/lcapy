"""This module provides continuous-time filter support.

Copyright 2022-2023 Michael Hayes, UCECE

"""

from .expr import expr, equation, ExprTuple
from .differentialequation import DifferentialEquation
from .functions import Derivative, Function, exp, cos
from .texpr import TimeDomainExpression
from .transfer import transfer
from .symbols import t, s, omega0, j, pi, f
from .utils import isiterable
import sympy as sym


class LTIFilter(object):

    def __init__(self, b, a):
        """Create continuous-time filter where `b` is a list or array of
        numerator coefficients and `a` is a list or array of
        denominator coefficients."""

        # For numerical stability it might be better to store filter
        # as ZPK form.

        if not isiterable(b):
            b = (b, )
        if not isiterable(a):
            a = (a, )

        self.a = ExprTuple(a)
        self.b = ExprTuple(b)

    def __repr__(self):
        """This is called by repr(expr).  It is used, e.g., when printing
        in the debugger."""

        return '%s(%s, %s)' % (self.__class__.__name__, self.b, self.a)

    @classmethod
    def from_ZPK(cls, Z, P, K=1):
        """Create LTIFilter given transfer function in zero-pole-gain form
        where `Z` is a list of zeroes, `P` is a list of poles, and `K`
        is a constant."""

        if not isiterable(Z):
            Z = (Z, )
        if not isiterable(P):
            P = (P, )

        N = expr(K)
        for z in Z:
            N *= (s - z)
        b = N.coeffs()

        D = expr(1)
        for p in P:
            D *= (s - p)
        a = D.coeffs()

        return cls(b, a)

    def transfer_function(self):
        """Return continuous-time impulse response (transfer function) in
        s-domain."""

        from .sym import ssym

        a = self.a.sympy
        b = self.b.sympy

        Nl = len(a)
        Nr = len(b)

        # Numerator of H(s)
        num = sym.S.Zero
        for i in range(Nr):
            num += b[i] * ssym**(Nr - i - 1)

        # Denominator for H(s)
        denom = sym.S.Zero
        for k in range(Nl):
            denom += a[k] * ssym**(Nl - k - 1)

        # Collect with respect to positive powers of the variable s
        num = sym.collect(sym.expand(num), ssym)
        denom = sym.collect(sym.expand(denom), ssym)

        Hs = transfer(sym.simplify(num / denom))
        return Hs

    def impulse_response(self):
        """Return continuous-time impulse response (transfer function) in time
        domain."""

        H = self.transfer_function()
        return H(t)

    def step_response(self):
        """Return continuous-time step response in time domain."""

        H = self.transfer_function()
        G = H / s

        return G(t)

    def frequency_response(self, **assumptions):
        """Return frequency response."""

        return self.transfer_function().frequency_response(**assumptions)

    def angular_frequency_response(self, **assumptions):
        """Return angular frequency response."""

        return self.transfer_function().angular_frequency_response(**assumptions)

    def group_delay(self):

        H = self.frequency_response()
        phi = H.phase
        tau = - Derivative(phi, f) / (2 * pi)
        return tau.simplify().general()

    def differential_equation(self, inputsym='x', outputsym='y'):
        """Return differential equation."""

        rhs = 0 * t

        for m, bn in enumerate(self.b):
            rhs += bn * Derivative(Function(inputsym)(t), t, m)

        for m, an in enumerate(self.a[1:]):
            rhs -= an * Derivative(Function(outputsym)(t), t, m + 1)

        lhs = self.a[0] * expr('y(t)')
        return DifferentialEquation(lhs, rhs, inputsym, outputsym)

    def sdomain_initial_response(self, ic=None):
        """Return s-domain response due to initial conditions.
           ic : list or tuple of initial conditions y(0), y'(0), y''(0), ...

        """

        from .sym import ssym

        if ic is None:
            ic = ()

        if not isiterable(ic):
            ic = (ic, )

        ic = ExprTuple(ic)

        Nl = len(self.a)
        if len(ic) != Nl - 1:
            raise ValueError('Expecting %d initial conditions, got %d'
                             % (Nl - 1, len(ic)))

        # Denominator for Yi(s)
        denom = self.a[0] * s**0
        num = 0 * s
        for k in range(1, Nl):
            aas = self.a[k] * s**(-k)
            denom += aas
            # Numerator for Yi(s)
            y0 = 0 * s
            for i in range(0, k):
                y0 += ic[i] * s**(i + 1)
            num += aas * y0

        # Collect with respect to positive powers of the variable s
        num = num.sympy * ssym**Nl
        denom = denom.sympy * ssym**Nl

        num = sym.collect(sym.expand(num), ssym)
        denom = sym.collect(sym.expand(denom), ssym)

        Ysi = expr(-sym.simplify(num / denom))
        Ysi.is_causal = True

        return Ysi

    def initial_response(self, ic=None):
        """Return response due to initial conditions in the time domain.
           ic : list of initial conditions y(0), y'(0), y''(0), ...
        """

        Ysi = self.sdomain_initial_response(ic)
        return Ysi(t)

    def response(self, x, ic=None, ni=None):
        """Calculate response of filter to input `x` given a list of initial conditions
        `ic` for time indexes specified by `ni`.  If `ni` is a tuple,
        this specifies the first and last (inclusive) time index.

        The initial conditions are valid prior to the time indices given by the ni
        `x` can be an expression, a sequence, or a list/array of values.
        """

        return 0

    def subs(self, *args, **kwargs):

        a = self.a.subs(*args, **kwargs)
        b = self.b.subs(*args, **kwargs)

        return self.__class__(b, a)

    def pdb(self):
        """Enter the python debugger."""

        import pdb
        pdb.set_trace()
        return self

    @property
    def is_stable(self):
        """Return True if all the poles of the filters transfer function
        are within the unit circle.

        See also is_marginally_stable."""

        return self.transfer_function().is_stable

    @property
    def is_marginally_stable(self):
        """Return True if all the poles of the filters transfer function
        are within or on the unit circle.

        See also is_stable."""

        return self.transfer_function().is_marginally_stable

    @classmethod
    def bessel(cls, order=2, kind='lpf', omega0=omega0):
        """Create a Bessel filter of specified order
        and angular cut-off frequency omega0."""

        from math import factorial

        D = expr(0)
        n = order
        for k in range(0, n + 1):

            a = factorial(2 * n - k) / (2 ** (n - k) *
                                        factorial(k) * factorial(n - k))
            D += a * (s / omega0)**k

        N = D.limit(s, 0)

        if kind == 'lpf':
            b = N.coeffs()
            a = D.coeffs()
        elif kind == 'hpf':
            H = (N / D).subs(s, omega0**2 / s)
            a = H.a
            b = H.b
        else:
            raise ValueError(
                'Bessel filter of kind %s not supported' % kind)

        return cls(b, a)

    @classmethod
    def butterworth(cls, order=2, kind='lpf', omega0=omega0):
        """Create a Butterworth filter of specified order
        and angular cut-off frequency omega0."""

        N = expr(1)
        D = expr(1)
        for k in range(1, order + 1):
            sk = omega0 * exp(j * (2 * k + order - 1) * pi / (2 * order))
            N *= omega0
            # Rewrite exp as cos to improve simplification.  For some reason
            # sympy does not simplify exp(-j pi / 4).
            sk = sk.rewrite(cos)
            D *= (s - sk)

        # Can convert LPF to HPF with s = omega0**2 / s

        if kind == 'lpf':
            b = N.coeffs()
        elif kind == 'hpf':
            b = [0] * (order + 1)
            b[0] = 1
        else:
            raise ValueError(
                'Butterworth filter of kind %s not supported' % kind)

        a = D.coeffs()

        return cls(b, a)


Bessel = LTIFilter.bessel
Butterworth = LTIFilter.butterworth

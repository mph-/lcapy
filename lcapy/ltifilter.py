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

    @classmethod
    def from_transfer_function(cls, H, normalize_a0=True):
        """Create LTIFilter given a transfer function."""

        H = transfer(H)
        if not H.is_rational_function:
            raise ValueError("Transfer function is not a rational function")

        N, D = H.as_N_D()

        bn = N.coeffs()
        an = D.coeffs()

        if normalize_a0:
            bn = [bx / an[0] for bx in bn]
            an = [ax / an[0] for ax in an]

        return cls(bn, an)

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

    def phase_response(self):

        H = self.frequency_response()
        return H.phase

    def group_delay(self):

        phi = self.phase_response()
        tau = -Derivative(phi, f) / (2 * pi)
        return tau.simplify().general()

    def phase_delay(self):

        phi = self.phase_response()
        tau = -phi / (2 * pi * f)
        return tau.simplify().general()

    def differential_equation(self, inputsym='x', outputsym='y'):
        """Return differential equation."""

        rhs = 0 * t
        lhs = 0 * t

        x = Function(inputsym)(t)
        y = Function(outputsym)(t)

        for m, bn in enumerate(reversed(self.b)):
            rhs += bn * Derivative(x, t, m)

        for m, an in enumerate(reversed(self.a)):
            lhs += an * Derivative(y, t, m)

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

    def discretize(self, method='bilinear', alpha=0.5):
        """Convert to a discrete-time linear time-invariant filter.

        The default method is 'bilinear'.  Other methods are:
        'impulse-invariance' 'bilinear', 'tustin', 'trapezoidal'
        'generalized-bilinear', 'gbf' controlled by the parameter
        `alpha`, 'euler', 'forward-diff', 'forward-euler'
        'backward-diff', 'backward-euler' 'simpson', 'matched-Z',
        'zero-pole-matching'"""

        H = self.transfer_function()

        return H.discretize(method, alpha).dlti_filter()

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
    def bessel(cls, N=2, Wn=omega0, btype='lowpass'):
        """Create a Bessel-Thomson (maximally flat delay) filter of specified
        order `N` and angular cut-off frequency Wn.

        """

        from math import factorial

        # TODO, add band-pass and band-stop with list for Wn.
        omega0 = Wn

        denom = expr(0)
        for k in range(0, N + 1):

            a = factorial(2 * N - k) / (2 ** (N - k) *
                                        factorial(k) * factorial(N - k))
            denom += a * (s / omega0)**k

        numer = denom.limit(s, 0)

        if btype == 'lowpass':
            b = numer.coeffs()
            a = denom.coeffs()
        elif btype == 'highpass':
            H = (numer / denom).subs(s, omega0**2 / s)
            a = H.a
            b = H.b
        else:
            raise ValueError(
                'Bessel filter of btype %s not supported' % btype)

        return cls(b, a)

    @classmethod
    def butterworth(cls, N=2, Wn=omega0, btype='lowpass'):
        """Create a Butterworth (maximally flat amplitude) filter of specified
        order `N` and angular cut-off frequency Wn.

        """

        # TODO, add band-pass and band-stop with list for Wn.
        omega0 = Wn

        numer = expr(1)
        denom = expr(1)
        for k in range(1, N + 1):
            sk = omega0 * exp(j * (2 * k + N - 1) * pi / (2 * N))
            numer *= omega0
            # Rewrite exp as cos to improve simplification.  For some reason
            # SymPy does not simplify exp(-j * pi / 4).
            sk = sk.rewrite(cos)
            denom *= (s - sk)

        # Can convert low-pass to high-pass with s = omega0**2 / s

        if btype == 'lowpass':
            b = numer.coeffs()
        elif btype == 'highpass':
            b = [0] * (N + 1)
            b[0] = 1
        else:
            raise ValueError(
                'Butterworth filter of btype %s not supported' % btype)

        a = denom.coeffs()

        return cls(b, a)


Bessel = LTIFilter.bessel
Butterworth = LTIFilter.butterworth

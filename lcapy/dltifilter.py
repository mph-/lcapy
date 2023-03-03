"""This module provides discrete-time filter support.

Copyright 2021--2023 Michael Hayes, UCECE

"""

from .expr import expr, equation, ExprTuple
from .functions import Function
from .nexpr import DiscreteTimeDomainExpression
from .differenceequation import DifferenceEquation
from .discretetime import n, z, seq
from .sequence import Sequence
from .utils import isiterable
from .sym import oo
import sympy as sym


class DLTIFilter(object):

    def __init__(self, b, a):
        """Create discrete-time filter where `b` is a list or array of
        numerator coefficients and `a` is a list or array of
        denominator coefficients.
        The order in a and b is according to z^0, z^(-1), ...

        Note, a recursive (IIR) filter with more than 3 coefficients
        for `a` is sensitive to truncation of the cofficients and can
        easily become unstable.  It is better to use a sequence
        of biquadratic filters.

        """

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

    def transfer_function(self):
        """Return discrete-time impulse response (transfer function) in
        z-domain."""

        from .sym import zsym

        a = self.a.sympy
        b = self.b.sympy

        Nl = len(a)
        Nr = len(b)

        # Numerator of H(z)
        num = sym.S.Zero
        for i in range(Nr):
            num += b[i] * zsym**(-i)

        # Denominator for H(z)
        denom = a[0] * zsym**0
        for k in range(1, Nl):
            az = a[k] * zsym**(-k)
            denom += az

        # Collect with respect to positive powers of the variable z
        num = sym.collect(sym.expand(num * zsym**Nl), zsym)
        denom = sym.collect(sym.expand(denom * zsym**Nl), zsym)

        Hz = expr(sym.simplify(num / denom))
        Hz.is_causal = True
        return Hz

    def impulse_response(self):
        """Return discrete-time impulse response (transfer function) in time domain."""

        H = self.transfer_function()
        return H(n)

    def step_response(self):
        """Return discrete-time step response in time domain."""

        H = self.transfer_function()
        G = (H * z) / (z - 1)
        G.is_causal = H.is_causal

        return G(n)

    def frequency_response(self, var=None, images=oo, **assumptions):
        """Return frequency response."""

        return self.transfer_function().DTFT(var, images, **assumptions)

    def difference_equation(self, inputsym='x', outputsym='y'):
        """
             difference_equation(inputsym='x', outputsym='y')
             Return difference equation.
             y[n] is a function of x[n-i] and y[n-i]     i>=0
        """

        rhs = 0 * n

        for m, bn in enumerate(self.b):
            rhs += bn * Function(inputsym)(n - m)

        for m, an in enumerate(self.a[1:]):
            rhs -= an * Function(outputsym)(n - (m + 1))

        lhs = self.a[0] * expr('y(n)')
        return DifferenceEquation(lhs, rhs, inputsym, outputsym)

    def zdomain_initial_response(self, ic=None, xic=None, left=True,
                                 inputsym='x', outputsym='y'):
        """
           Return z-domain response due to initial conditions.
           zdomain_initial_response(ic=None, xic=None, left=True, inputsym='x', outputsym='y')
           left=True (default)
              ic : list with initial values y[-1], y[-2], ...,
              xic : list with initial values x[-1], x[-2], ...,
           left=False
              ic : list with initial values y[0], y[1], ...
              xic : list with initial values x[0], x[1], ...
           If no ic or xic are given symbolic expressions are used
        """

        from .sym import zsym

        Nl = len(self.a)

        # default values for ic
        if ic is None:
            ic = []
            n_pos = range(-1, -Nl, -1)
            if not left:
                n_pos = range(0, Nl - 1)
            for k in n_pos:
                ic += [sym.Symbol('%s[%i]' % (outputsym, k))]
            ic = tuple(ic)

        # default values for xic
        if xic is None:
            xic = []
            n_pos = range(-1, -Nl, -1)
            if not left:
                n_pos = range(0, Nl - 1)
            for k in n_pos:
                xic += [sym.Symbol('%s[%i]' % (inputsym, k))]
            xic = tuple(xic)

        if not isiterable(ic):
            ic = (ic, )
        if not isiterable(xic):
            xic = (xic, )

        Nl = len(self.a)
        if len(ic) != Nl - 1:
            raise ValueError('Expecting %d initial conditions for ic, got %d'
                             % (Nl - 1, len(ic)))
        if len(xic) != Nl - 1:
            raise ValueError('Expecting %d initial conditions for xic, got %d'
                             % (Nl - 1, len(ic)))

        # Initial condition y[-1], y[-2], .....
        if left:
            # Denominator for Yi(z)
            denom = self.a[0] * z**0
            num = 0 * z
            for k in range(1, Nl):
                az = self.a[k] * z**(-k)
                denom += az
                # Numerator for Yi(z)
                y0 = 0 * z
                x0 = 0 * z
                for i in range(0, k):
                    y0 += ic[i] * z**(i + 1)
                    x0 += xic[i] * z**(i + 1)
                # ic conditions
                num -= az * y0
                # xic conditions
                if len(self.b) >= k + 1:
                    num += x0 * self.b[k] * z**(-k)

            # Collect with respect to positive powers of the variable z
            num = num.sympy * zsym**(Nl - 1)
            denom = denom.sympy * zsym**(Nl - 1)

        # Initial condition y[0], y[1], .....
        else:
            # Denominator for Yi(z)
            denom = self.a[Nl - 1]
            num = 0 * z
            for k in range(0, Nl - 1):
                az = self.a[k] * z**(Nl - 1 - k)
                denom += az
                # Numerator for Yi(z)
                y0 = 0 * z
                x0 = 0 * z
                for i in range(0, Nl - 1 - k):
                    y0 += ic[i] * z ** (-i)
                    x0 -= xic[i] * z ** (-i)
                # ic conditions
                num += az * y0
                # xic conditions
                if len(self.b) >= k + 1:
                    num += x0 * self.b[k] * z**(Nl - 1 - k)

        num = sym.collect(sym.expand(num), zsym)
        denom = sym.collect(sym.expand(denom), zsym)

        Yzi = expr(num / denom)
        Yzi.is_causal = True

        return Yzi

    def initial_response(self, ic=None, xic=None, left=True, inputsym='x', outputsym='y'):
        """
           initial_response(ic=None, xic=None, left=True, inputsym='x', outputsym='y')
           Return response due to initial conditions in the time domain.
           left=True (default)
              ic : list with initial values y[-1], y[-2], ...
              xic : list with initial values x[-1], x[-2], ...
           left=False
              ic : list with initial values y[0], y[1], ...
              xic : list with initial values x[0], x[1], ...
           If no ic or xic are given symbolic expressions are used
        """

        Yzi = self.zdomain_initial_response(
            ic=ic, xic=xic, left=left, inputsym=inputsym, outputsym=outputsym)
        return Yzi(n)

    def response(self, x, ic=None, ni=None):
        """
        response(x, ic=None, ni=None)
        Calculate response of filter to input `x` at time indexes specified by `ni`,
        given a list of initial conditions `ic`.

        - `x` can be an expression, a sequence, or a list/array of values.
        - If `ni` is a tuple, this specifies the first and last (inclusive) time index.
          Otherwise, it should be list of integers, monotonically increasing by 1.
          If `ni` is `None`, or unspecified, it defaults to `(-len(ic), 10)`.
        - `ic` specifies the initial conditions, `[y[-1], y[-2], ...]`.
        - If `ic` is `None` the initial conditions are assumed to be zero.
        """

        from numpy import arange, ndarray

        if (ic is None and ni is None and isinstance(x, ndarray)
                and self.a.symbols == {} and self.b.symbols == {}):

            from scipy.signal import lfilter

            # Calculate numerical response
            b = self.b.fval
            a = self.a.fval
            return lfilter(b, a, x)

        if ic is None:
            ic = [0] * (len(self.a) - 1)

        if not isiterable(x) and not isinstance(x, DiscreteTimeDomainExpression):
            x = (x, )

        if not isiterable(ic):
            ic = (ic, )

        if isinstance(x, (tuple, list, ndarray)):
            x = seq(x)
        elif not isinstance(x, (Sequence, DiscreteTimeDomainExpression)):
            raise ValueError(
                'The input x must be a scalar, tuple, sequence, nexpr, list, or array')

        Ni = len(ic)

        if Ni != len(self.a) - 1:
            raise ValueError(
                "Expected %d initial conditions, got %d" % (len(self.a) - 1, Ni))

        if ni is None:
            ni = (-Ni, 10)

        if isinstance(ni, tuple):
            ni = arange(ni[0], ni[1] + 1)

        Nn = ni[-1] + 1

        # Order right hand side
        Nr = len(self.b)

        y_tot = list(ic[-1::-1]) + [0] * Nn

        a_r = self.a[-1:-1 - Ni:-1]
        for i in range(Nn):
            # Get previous y vals (sliding window)
            pre_y = y_tot[i:i + Ni]

            # Calculate rhs of new value
            if isinstance(x, Sequence):
                rhs = sum(self.b[l] * x[i - l] for l in range(Nr))
            else:
                rhs = sum(self.b[l] * x(i - l) for l in range(Nr))

            # Add lhs
            y_tot[i + Ni] = -1 / self.a[0] * \
                sum(csi * ysi for csi, ysi in zip(a_r, pre_y)) + \
                rhs / self.a[0]

        # Number of zeros to prepend
        Nz = 0
        if ni[0] < -Ni:
            Nz = -(Ni + ni[0])
            y_tot = Nz * [0] + y_tot

        ret_seq = seq(y_tot[ni[0] + Ni + Nz:], tuple(ni))

        return ret_seq

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

    def inverse(self):
        """Create inverse filter.  Note, in practice, an inverse filter
        is sensitive to noise.  It is better to use a Wiener filter."""

        return self.__class__(self.a, self.b)

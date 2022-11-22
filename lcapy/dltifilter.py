"""This module provides discrete-time filter support.

Copyright 2021--2022 Michael Hayes, UCECE

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
        """Return difference equation."""

        rhs = 0 * n

        for m, bn in enumerate(self.b):
            rhs += bn * Function(inputsym)(n - m)

        for m, an in enumerate(self.a[1:]):
            rhs -= an * Function(outputsym)(n - (m + 1))

        lhs = self.a[0] * expr('y(n)')
        return DifferenceEquation(lhs, rhs, inputsym, outputsym)

    def zdomain_initial_response(self, ic=None):
        """Return z-domain response due to initial conditions.
           ic : list with initial values y[-1], y[-2], ...
        """

        from .sym import zsym

        if ic is None:
            ic = ()

        if not isiterable(ic):
            ic = (ic, )

        Nl = len(self.a)
        if len(ic) != Nl:
            raise ValueError('Expecting %d initial conditions, got %d'
                             % (Nl, len(ic)))

        # Denominator for Yi(z)
        denom = self.a[0] * z**0
        num = 0 * z
        for k in range(1, Nl):
            az = self.a[k] * z**(-k)
            denom += az
            # Numerator for Yi(z)
            y0 = 0 * z
            for i in range(0, k):
                y0 += ic[i] * z**(i + 1)
            num += az * y0

        # Collect with respect to positive powers of the variable z
        num = num.sympy * zsym**Nl
        denom = denom.sympy * zsym**Nl

        num = sym.collect(sym.expand(num), zsym)
        denom = sym.collect(sym.expand(denom), zsym)

        Yzi = expr(-sym.simplify(num / denom))
        Yzi.is_causal = True

        return Yzi

    def initial_response(self, ic=None):
        """Return response due to initial conditions in the time domain.
           ic : list with initial values y[-1], y[-2], ...
        """

        Yzi = self.zdomain_initial_response(ic)
        return Yzi(n)

    def response(self, x, ic=None, ni=None):
        """Calculate response of filter to input `x` given a list of initial conditions
        `ic` for time indexes specified by `ni`.  If `ni` is a tuple,
        this specifies the first and last (inclusive) time index.

        The initial conditions are valid prior to the time indices given by the ni
        `x` can be an expression, a sequence, or a list/array of values.
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

        NO = len(ic)

        if NO != len(self.a) - 1:
            raise ValueError(
                "Expected %d initial conditions, got %d" % (len(self.a) - 1, NO))

        if ni is None:
            ni = (0, 10)

        if isinstance(ni, tuple):
            ni = arange(ni[0], ni[1] + 1)

        Nn = len(ni)

        # Order right hand side
        Nr = len(self.b)

        y_tot = list(ic[-1::-1]) + Nn * [0]

        a_r = self.a[-1:-1-NO:-1]
        for i, nval in enumerate(ni):
            # Get previous y vals (sliding window)
            pre_y = y_tot[i:i + NO]

            # Calculate rhs of new value
            if isinstance(x, Sequence):
                rhs = sum(self.b[l] * x[nval - l] for l in range(Nr))
            else:
                rhs = sum(self.b[l] * x(nval - l) for l in range(Nr))

            # Add lhs
            y_tot[i + NO] = -1 / self.a[0] * \
                sum(csi * ysi for csi, ysi in zip(a_r, pre_y)) + \
                rhs / self.a[0]

        # Solution, without initial values
        ret_seq = seq(y_tot[NO:], ni)

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

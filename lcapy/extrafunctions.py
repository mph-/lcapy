"""This module defines additional functions to those defined by SymPy.
They are for internal use by Lcapy.  The user versions are defined in
function.py.

Copyright 2020--2022 Michael Hayes, UCECE

"""

import sympy as sym
from sympy.core import S, Integer
from sympy.core.logic import fuzzy_not
from .config import unitstep_zero


class Degrees(sym.Function):
    """Convert radians to degrees."""

    @classmethod
    def eval(cls, nval):
        """
        Convert radians to degrees.
        """

        return nval / sym.pi * 180


class Radians(sym.Function):
    """Convert degrees to radians."""

    @classmethod
    def eval(cls, nval):
        """
        Convert degrees to radians.
        """

        return nval / 180 * sym.pi


class UnitImpulse(sym.Function):
    """Discrete unit impulse function."""

    is_integer = True

    @classmethod
    def eval(cls, nval):
        """
        Evaluates the discrete unit impulse function.
        """

        if nval.is_zero:
            return S.One
        elif fuzzy_not(nval.is_zero):
            return S.Zero


class UnitStep(sym.Function):
    """Discrete unit step function."""

    is_integer = True

    @classmethod
    def eval(cls, nval, zero=None):
        """
        Evaluates the discrete unit step function.   This is defined
        as 1 for n >= 0 and 0 otherwise.
        """

        if nval.is_zero:
            if zero is None:
                zero = unitstep_zero
            return zero

        if nval.is_nonnegative:
            return S.One
        elif nval.is_negative:
            return S.Zero


class rect(sym.Function):
    """The rectangle function.

    rect(t) = 1 for abs(t) <= 0.5 otherwise 0."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the rectangle function.
        """

        if val.is_Number:
            if val < -0.5 or val > 0.5:
                return S.Zero
            return S.One

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return sym.Heaviside(x + S.Half) - sym.Heaviside(x - S.Half)


class dtrect(sym.Function):
    """The discrete-time rectangle function."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the rectangle function for discrete-time signals.
        """

        if val.is_Number:
            if val < -0.5 or val >= 0.5:
                return S.Zero
            return S.One

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return UnitStep(x + S.Half) - UnitStep(x - S.Half)


class dtsign(sym.Function):
    """The discrete-time signum function."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the signum function for discrete-time signals.
        """

        if val.is_Number:
            if val >= 0:
                return S.One
            return S.Zero


class sincn(sym.Function):
    """Normalized sinc function :math:`\sin(\pi x)/(\pi x)`."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the normalized sinc (cardinal sine) function.
        This is what NumPy uses but not SymPy.
        """

        if val.is_Number:
            if val == 0:
                return S.One
            x = sym.pi * val
            return sym.sin(x) / x

    def rewrite(self, *args, **hints):

        x = sym.pi * self.args[0]
        return sym.sin(x) / x


class sincu(sym.Function):
    """Unnormalized sinc function :math:`\sin(x)/(x)`."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the unnormalized sinc (cardinal sine) function.
        This is what SymPy uses but not NumPy.
        """

        if val.is_Number:
            if val == 0:
                return S.One
            return sym.sin(val) / val

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return sym.sin(x) / x


class psinc(sym.Function):
    """Periodic sinc function :math:`\sin(M * pi * x)/(M * sin(pi * x))`."""

    @classmethod
    def eval(cls, M, val):
        """
        Evaluates the periodic sinc function.
        """
        if val.is_Number:
            if val == 0:
                return S.One
            if val.is_integer and val.is_even:
                return S.One
            if val.is_integer and val.is_odd:
                return -S.One

            x = sym.pi * val
            return sym.sin(M * x) / (M * sym.sin(x))

    def rewrite(self, *args, **hints):

        M = self.args[0]
        x = sym.pi * self.args[1]
        # If evaluate, SymPy will convert sin(n * pi) to 0.
        return sym.sin(M * x, evaluate=False) / (M * sym.sin(x, evaluate=False))


class tri(sym.Function):
    """The triangle function.

    tri(t) = 1 - abs(t) for abs(t) <= 1 otherwise 0."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the triangle function.
        """

        if val.is_Number:

            if val >= 1:
                return S.Zero
            elif val <= -1:
                return S.Zero
            else:
                return 1 - abs(val)

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return ramp(x + 1) - 2 * ramp(x) + ramp(x - 1)


class trap(sym.Function):
    """The trapezoidal function."""

    @classmethod
    def eval(cls, val, alpha):
        """
        Evaluates the trapezoid function.   This is rect(t / alpha).convolve(rect(t)).
        trap(t, 0) = rect(t) and trap(t, 1) = tri(t).
        """

        if val.is_Number and alpha.is_number:
            absval = abs(val)
            foo = absval - 0.5

            if alpha == S.Zero:
                if val < -0.5 or val > 0.5:
                    return S.Zero
                return S.One

            if foo >= 0.5 * alpha:
                return S.Zero
            elif foo <= -0.5 * alpha:
                return S.One
            else:
                return 0.5 - foo / alpha

    def rewrite(self, *args, **hints):

        x = self.args[0]
        alpha = self.args[1]
        if alpha == 0:
            return rect(x)
        elif alpha == 1:
            return tri(x)
        return rampstep((x + (1 + alpha) / 2) / alpha) - rampstep((x - (1 - alpha) / 2) / alpha)


class ramp(sym.Function):
    """The ramp function.

    ramp(t) = t for t > 0 otherwise 0."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the ramp function.
        """

        if val.is_Number:

            if val >= 0:
                return val
            else:
                return 0

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return x * sym.Heaviside(x)


class rampstep(sym.Function):
    """The rampstep function.

    rampstep(t) = t for 0 < t < 1, 1 for t > 1 and otherwise 0."""

    @classmethod
    def eval(cls, val):
        """
        Evaluates the rampstep function.
        """

        if val.is_Number:

            if val >= 0 and val < 1:
                return val
            elif val > 1:
                return 1
            else:
                return 0

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return ramp(x) - ramp(x - 1)

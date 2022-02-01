"""This module provides support for inverse discrete-time Fourier
transforms.

 It calculates the inverse discrete-time Fourier transform using:

   x(n) = dt * \int_{-1 / (2 * dt)}^{1 / (2 * dt))} X(f) e^{j * 2 * \pi * n * dt * f} df

where dt is the sampling period.  The integral just needs to be
performed over any interval of length 1 / dt.

Copyright 2021--2022 Michael Hayes, UCECE

"""

import sympy as sym
from sympy import oo, DiracDelta
from .transformer import BilateralInverseTransformer
from .sym import sympify, AppliedUndef, j, pi
from .sym import dt
from .extrafunctions import UnitImpulse, UnitStep
from .utils import factor_const, scale_shift, remove_images
from .matrix import Matrix

__all__ = ('IDTFT', 'inverse_discrete_time_fourier_transform')


class IDTFTTransformer(BilateralInverseTransformer):

    name = 'inverse DTFT'

    def key(self, expr, f, n, **assumptions):
        return expr, f, n

    def noevaluate(self, expr, f, n):

        foo = expr * sym.exp(2 * j * pi * n * dt * f)
        result = dt * sym.Integral(foo, (f, -1 / (2 * dt), 1 / (2 * dt)))
        return result

    def check(self, expr, f, n, **assumptions):

        if expr.has(n):
            self.error('Expression depends on n')

    def rewrite(self, expr, var):
        """Remove images."""

        return remove_images(expr, var, dt)

    def sympy(self, expr, f, n):

        foo = expr * sym.exp(2 * j * pi * n * dt * f)
        result = dt * sym.integrate(foo, (f, -1 / (2 * dt), 1 / (2 * dt)))
        return result

    def func(self, expr, f, n):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], f)

        # Convert v(n) to V(f), etc.
        name = expr.func.__name__
        func = sym.Function(name[0].lower() + name[1:])

        result = func(n / scale) / abs(scale)

        if shift != 0:
            result = result * sym.exp(-2 * sym.I * sym.pi * f * shift / scale)

        return result

    def function(self, expr, f, n):

        # Handle expressions with a function of FOO, e.g.,
        # v(n), v(n) * y(n),  3 * v(n) / n, v(4 * a * n), etc.,

        if not expr.has(AppliedUndef):
            self.error()

        const, expr = factor_const(expr, f)

        if isinstance(expr, AppliedUndef):
            return self.func(expr, f, n) * const

        nsym = sympify(str(n))
        expr = expr.subs(nsym, n)

        rest = sym.S.One
        undefs = []
        for factor in expr.as_ordered_factors():
            if isinstance(factor, AppliedUndef):
                if factor.args[0] != f:
                    self.error('Weird function %s not of %s' % (factor, f))
                undefs.append(factor)
            else:
                rest *= factor

        if rest.has(AppliedUndef):
            # Have something like 1/V(f)
            self.error()

        exprs = undefs
        if rest.has(n):
            exprs = exprs + [rest]
            rest = sym.S.One

        result = self.term(exprs[0], f, n) * rest

        if len(exprs) == 1:
            return result * const

        self.error('TODO')

    def term(self, expr, f, n):

        const, expr = factor_const(expr, f)

        # Check for constant.
        if not expr.has(f):
            return expr * UnitImpulse(n) * const

        if expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            return self.function(expr, f, n) * const

        # This handles DiracDelta the long way but could head off at pass...
        return const * self.sympy(expr, f, n)


idtft_transformer = IDTFTTransformer()


def inverse_discrete_time_fourier_transform(expr, f, n, evaluate=True,
                                            **assumptions):
    """Compute bilateral inverse discrete-time Fourier transform of expr.

    Undefined functions such as X(f) are converted to x(n).

    """

    return idtft_transformer.transform(expr, f, n, evaluate=evaluate,
                                       **assumptions)


def IDTFT(expr, f, n, evaluate=True, **assumptions):
    """Compute bilateral inverse discrete-time Fourier transform of expr.

    Undefined functions such as X(f) are converted to x(n).
    """

    return idtft_transformer.transform(expr, f, n, evaluate=evaluate,
                                        **assumptions)

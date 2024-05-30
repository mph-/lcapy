"""This module provides support for Hilbert transforms.  It calculates
the Hilbert transform using:


Copyright 2024 Michael Hayes, UCECE

"""

from sympy import pi, oo, exp, sin, cos, log, sign, Integral, I, DiracDelta
from .sym import symsimplify, j, tausym
from .transformer import BilateralForwardTransformer
from .utils import factor_const
from .extrafunctions import rect

__all__ = ('HT', )


# If g(t) has a Hilbert transform H g(t) then
# g(t - t0) has a Hilbert transform H g(t - t0) and
# g(a t)has a Hilbert transform sign(a) H g(a t).
#
# H(g(t) * h(t)) = H g(t) * h(t) = g(t) * H h(t).
#
# H dg(t)/dt = d H g(t) dt




class HilbertTransformer(BilateralForwardTransformer):

    name = 'Hilbert transform'

    def key(self, expr, t, f, **kwargs):
        return expr, t, f

    def simplify_term(self, expr, var):

        return symsimplify(expr)

    def term(self, expr, t, f=None):

        const, expr = factor_const(expr, t)

        if self.is_inverse:
            const = -const

        # Check for constant.
        if not expr.has(t):
            return 0

        if expr.is_Pow and expr.args[1] == -1:
            # Handle 1 / t
            arg = expr.args[0]
            if arg.is_polynomial(t) and arg.as_poly(t).is_linear:
                c1 = arg.expand().coeff(t, 1)

                return pi * DiracDelta(arg) / c1

        if len(expr.args) == 0:
            return const * self.noevaluate(expr, t)

        arg = expr.args[0]

        if not arg.is_polynomial(t) or not arg.as_poly(t).is_linear:
            return const * self.noevaluate(expr, t)

        if not expr.is_Function:
            return const * self.noevaluate(expr, t)

        c1 = arg.expand().coeff(t, 1)

        if expr.func == DiracDelta:
            return const * sign(c1) * 1 / (pi * arg)

        if expr.func == cos:
            return const * sign(c1) * sin(arg)

        elif expr.func == sin:
            return const * sign(c1) * -cos(arg)

        elif expr.func == rect:
            return const * sign(c1) * -log(abs((arg - 1 / 2) / (arg + 1 / 2))) / pi

        elif expr.func == exp and c1.has(I):
            return const * sign(c1 / I) * -I * expr

        # TODO, add sinc pulse, Cauchy pulse, 1 / t


        q =  Integral(expr / (pi * (t - tausym)), (tausym, -oo, oo))
        return const * q.principal_value()

    def noevaluate(self, expr, t, f=None):

        return Integral(expr / (pi * (t - tausym)), (tausym, -oo, oo))

    def check(self, expr, t, f=None):
        pass

    def rewrite(self, expr, var):
        return expr


hilbert_transformer = HilbertTransformer()


def HT(expr, t, **kwargs):
    """Compute Hilbert transform of expr."""

    return hilbert_transformer.transform(expr, t, **kwargs)


def hilbert_transform(expr, t, **kwargs):
    """Compute Hilbert transform of expr."""

    return hilbert_transformer.transform(expr, t, **kwargs)

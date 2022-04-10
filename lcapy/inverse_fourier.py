"""This module provides support for inverse Fourier transforms.  It
calculates the bilateral inverse Fourier transform using:

   s(t) = \int_{-\infty}^{\infty} S(f) e^{j * 2 * \pi * t} df

It also allows functions that strictly do not have an inverse Fourier
transform by using Dirac deltas.

Copyright 2016--2021 Michael Hayes, UCECE
"""

from .fourier import FourierTransformer
from .sym import j, pi
import sympy as sym


class InverseFourierTransformer(FourierTransformer):

    name = 'inverse Fourier transform'
    is_inverse = True

    def noevaluate(self, expr, f, t):

        return sym.Integral(expr * sym.exp(j * 2 * pi * f * t), (f, -sym.oo, sym.oo))

    def check(self, expr, f, t):

        if expr.has(t):
            self.error('Expression depends on t')


inverse_fourier_transformer = InverseFourierTransformer()


def IFT(expr, f, t, **kwargs):
    """Compute bilateral inverse Fourier transform of expr.

    Undefined functions such as V(f) are converted to v(t)

    This also handles some expressions that do not really have an
    inverse Fourier transform, such as a, cos(a * f), sin(a * f), exp(I *
    a * f)."""

    return inverse_fourier_transformer.transform(expr, f, t, **kwargs)


def inverse_fourier_transform(expr, f, t, **kwargs):
    """Compute bilateral inverse Fourier transform of expr.

    Undefined functions such as V(f) are converted to v(t)

    This also handles some expressions that do not really have an
    inverse Fourier transform, such as a, cos(a * f), sin(a * f), exp(I *
    a * f)."""

    return inverse_fourier_transformer.transform(expr, f, t, **kwargs)

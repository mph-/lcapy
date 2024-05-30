"""This module provides support for inverse Hilbert transforms.

Copyright 2024 Michael Hayes, UCECE
"""

from .hilbert import HilbertTransformer
from .sym import j, pi
import sympy as sym


class InverseHilbertTransformer(HilbertTransformer):

    name = 'inverse Hilbert transform'
    is_inverse = True

    def noevaluate(self, expr, f, t):

        return sym.Integral(expr * sym.exp(j * 2 * pi * f * t), (f, -sym.oo, sym.oo))


inverse_hilbert_transformer = InverseHilbertTransformer()


def IHT(expr, t, **kwargs):
    """Compute inverse Hilbert transform of expr.

    Note, if the original function has a constant then
    it cannot be recovered."""

    return inverse_hilbert_transformer.transform(expr, t, **kwargs)


def inverse_hilbert_transform(expr, t, **kwargs):
    """Compute inverse Hilbert transform of expr.

    Note, if the original function has a constant then
    it cannot be recovered."""

    return inverse_hilbert_transformer.transform(expr, t, **kwargs)

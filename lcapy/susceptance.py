"""This module provides susceptance support.

Copyright 2021--2022 Michael Hayes, UCECE

"""
from .expr import expr
from .sym import j
from .units import u as uu
from warnings import warn


def susceptance(arg, **assumptions):
    """Create an admittance class for the specified susceptance.

    The susceptance, B, is defined as the imaginary part of the admittance.

    Y(omega) = G(omega) + j * B(omega)"""

    expr1 = expr(arg, frequency=True, **assumptions)
    if expr1.is_laplace_domain:
        warn('Specifying Laplace domain for susceptance: %s' % expr1)

    if expr1.is_imaginary:
        warn('Susceptance %s should be real' % expr1)

    expr1.units = uu.siemens
    return expr1

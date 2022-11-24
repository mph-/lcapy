"""This module provides conductance support.

Copyright 2021--2022 Michael Hayes, UCECE

"""
from .expr import expr
from .units import u as uu
from warnings import warn


def conductance(arg, **assumptions):
    """Create an admittance class for the specified conductance.

    The conductance, G, is defined as the real part of the admittance:

    Y(omega) = G(omega) + j B(omega)"""

    expr1 = expr(arg, frequency=True, **assumptions)
    if expr1.is_imaginary:
        warn('Conductance %s should be real' % expr1)

    expr1.units = uu.siemens
    return expr1

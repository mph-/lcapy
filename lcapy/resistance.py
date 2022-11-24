"""This module provides resistance support.

Copyright 2021--2022 Michael Hayes, UCECE

"""
from .expr import expr
from .units import u as uu
from warnings import warn


def resistance(arg, **assumptions):
    """Create an impedance class for the specified resistance.

    The resistance, R, is defined as the real part of the impedance:

    Z(omega) = R(omega) + j * X(omega)"""

    expr1 = expr(arg, frequency=True, **assumptions)
    if expr1.is_imaginary:
        warn('Resistance %s should be real' % expr1)

    expr1.units = uu.ohms
    return expr1

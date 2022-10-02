"""This module provides reactance support.

Copyright 2021--2022 Michael Hayes, UCECE

"""
from .expr import expr
from .sym import j
from warnings import warn


def reactance(arg, **assumptions):
    """Create an impedance class for the specified reactance.

    The reactance, X, is defined as the imaginary part of the impedance:

    Z(omega) = R(omega) + j * X(omega)"""

    expr1 = expr(arg, **assumptions)
    if expr1.is_laplace_domain:
        warn('Specifying Laplace domain for reactance: %s' % expr1)

    if expr1.is_imaginary:
        warn('Reactance %s should be real' % expr1)

    expr1 = j * expr1

    try:
        expr1 = expr1.as_impedance()
    except:
        raise ValueError('Cannot represent %s(%s) as impedance' %
                         (expr1.__class__.__name__, expr1))

    return expr1

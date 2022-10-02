"""This module provides conductance support.

Copyright 2021--2022 Michael Hayes, UCECE

"""
from .expr import expr
from warnings import warn


def conductance(arg, **assumptions):
    """Create an admittance class for the specified conductance.

    The conductance, G, is defined as the real part of the admittance:

    Y(omega) = G(omega) + j B(omega)"""

    expr1 = expr(arg, **assumptions)
    if expr1.is_imaginary:
        warn('Conductance %s should be real' % expr1)

    try:
        expr1 = expr1.as_admittance()
    except:
        raise ValueError('Cannot represent %s(%s) as admittance' %
                         (expr1.__class__.__name__, expr1))

    return expr1

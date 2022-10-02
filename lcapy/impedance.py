"""This module provides impedance support.

Copyright 2019--2022 Michael Hayes, UCECE

"""
from __future__ import division
from .expr import expr
from .sexpr import LaplaceDomainExpression


def impedance(arg, causal=True, **assumptions):
    """Create an impedance class for the specified impedance.

    Z(omega) = R(omega) + j * X(omega)

    where G is the conductance and B is the susceptance.

    Impedance is the reciprocal of admittance:

    Z(omega) = 1 / Y(omega)

    """

    expr1 = expr(arg, causal=causal, **assumptions)

    try:
        expr1 = expr1.as_impedance()
    except:
        raise ValueError('Cannot represent %s(%s) as impedance' %
                         (expr1.__class__.__name__, expr1))

    return expr1

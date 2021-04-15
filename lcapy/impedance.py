"""This module provides impedance support.

Copyright 2019--2020 Michael Hayes, UCECE

"""
from __future__ import division
from .expr import expr
from .sexpr import LaplaceDomainExpression


def impedance(arg, causal=True, **assumptions):
    """Generic impedance factory function.

    Y(omega) = G(omega) + j * B(omega)

    where G is the conductance and B is the susceptance.

    Admittance is the reciprocal of impedance,

    Z(omega) = 1 / Y(omega)

    """

    expr1 = expr(arg, causal=causal, **assumptions)    

    try:
        expr1 = expr1.as_impedance()
    except:    
        raise ValueError('Cannot represent %s(%s) as impedance' % (expr1.__class__.__name__, expr1))
        
    return expr1
    

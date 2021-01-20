"""This module provides admittance support.

Copyright 2019-2020 Michael Hayes, UCECE

"""
from __future__ import division
from .expr import expr
from .sexpr import LaplaceDomainExpression

    
def admittance(arg, **assumptions):
    """Generic admittance factory function.

    Y(omega) = G(omega) + j * B(omega)

    where G is the conductance and B is the susceptance.

    Admittance is the reciprocal of impedance,

    Z(omega) = 1 / Y(omega)

    """    

    expr1 = expr(arg, **assumptions)

    try:
        expr1 = expr1.as_admittance()
    except:
        raise ValueError('Cannot represent %s(%s) as admittance' % (expr1.__class__.__name__, expr1))        

    return expr1


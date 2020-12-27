"""This module provides transfer function support.

Copyright 2020 Michael Hayes, UCECE

"""
from .expr import expr
from .units import u as uu

def transfer(arg, **assumptions):

    expr1 = expr(arg, **assumptions)
    if expr1.is_admittance:
        return expr1.apply_unit(1 / uu.ohms)

    try:
        expr1 = expr1.as_transfer()
    except:
        raise ValueError('Cannot represent %s(%s) as transfer function' % (expr1.__class__.__name__, expr1))

    return expr1


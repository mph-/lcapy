"""This module provides transfer function support.

Copyright 2020--2022 Michael Hayes, UCECE

"""
from .expr import expr
from .units import u as uu


def transfer(numer, denom=1, causal=True, **assumptions):
    """Create transfer function assuming zero initial conditions
    from the ratio `numer / denom`."""

    numer = expr(numer, causal=causal, **assumptions)
    denom = expr(denom, causal=causal, **assumptions)

    numer = numer.zero_initial_conditions()
    denom = denom.zero_initial_conditions()
    expr1 = numer / denom

    if expr1.is_admittance:
        return expr1.apply_unit(1 / uu.ohms)

    try:
        expr1 = expr1.as_transfer()
    except:
        raise ValueError('Cannot represent %s(%s) as transfer function' % (
            expr1.__class__.__name__, expr1))

    return expr1

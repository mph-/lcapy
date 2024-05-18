"""This module provides transfer function support.

Copyright 2020--2024 Michael Hayes, UCECE

"""
from .expr import expr
from .units import u as uu
from .functions import Mul, Pow


def transfer(numer, denom=1, causal=True, **assumptions):
    """Create transfer function assuming zero initial conditions
    from the ratio `numer / denom`."""

    numer = expr(numer, causal=causal, **assumptions)
    denom = expr(denom, causal=causal, **assumptions)

    numer = numer.zero_initial_conditions()
    denom = denom.zero_initial_conditions()

    if False:
        expr1 = numer / denom
    else:
        # Avoid pole-zero cancellation by deferring evaluation
        expr1 = Mul(numer, Pow(denom, -1), evaluate=False)

    if expr1.is_admittance:
        expr1.units = 1 / uu.ohms
        return expr1

    try:
        expr1 = expr1.as_transfer()
    except:
        raise ValueError('Cannot represent %s(%s) as transfer function' % (
            expr1.__class__.__name__, expr1))

    return expr1

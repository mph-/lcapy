"""This module contains functions for approximating expressions.

Copyright 2021 Michael Hayes, UCECE

"""

from .simplify import expand_hyperbolic_trig
from sympy import exp


def approximate_fractional_power(expr, method='pade', order=2):
    """This is an experimental method to approximate s**a, where a is
    fractional, with a rational function using a Pade approximant.

    """

    if method != 'pade':
        raise ValueError('Method %s unsupported, must be pade')

    v = expr.var

    def query(expr):

        if not expr.is_Pow:
            return False
        if expr.args[0] != v:
            return False
        if expr.args[1].is_Number and not expr.args[1].is_Integer:
            return True
        if expr.args[1].is_Symbol and not expr.args[1].is_Integer:
            return True
        return False

    def value1(expr):

        a = expr.args[1]
        n = v * (a + 1) + (1 - a)
        d = v * (a - 1) + (1 + a)
        return n / d

    def value2(expr):

        a = expr.args[1]
        n = v**2 * (a**2 + 3 * a + 2) + v * (8 - a**2) + (a**2 - 3 * a + 2)
        d = v**2 * (a**2 - 3 * a + 2) + v * (8 - a**2) + (a**2 + 3 * a + 2)
        return n / d

    if order == 1:
        value = value1
    elif order == 2:
        value = value2
    else:
        raise ValueError('Can only handle order 1 and 2 at the moment')

    expr = expr.expr
    expr = expr.replace(query, value)

    return expr


def _approximate_exp_pade(expr, order=1):

    def query(expr):
        return expr.is_Function and expr.func == exp

    def value(expr):
        arg = expr.args[0]

        if order == 1:
            return (2 + arg) / (2 - arg)
        elif order == 2:
            return (12 + 6 * arg + arg**2) / (12 - 6 * arg + arg**2)

        from math import factorial

        numer = 0
        denom = 0
        for k in range(order + 1):
            scale = factorial(2 * order - k) // (factorial(k) * factorial(order - k))
            numer += scale * arg**k
            denom += scale * (-arg)**k
        return numer / denom

    return expr.replace(query, value)


def approximate_exp(expr, method='pade', order=1):
    """Approximate exp(a)."""

    if method != 'pade':
        raise ValueError('Method %s unsupported, must be pade')

    return _approximate_exp_pade(expr, order)


def approximate_hyperbolic_trig(expr, method='pade', order=1):
    """Approximate cosh(a), sinh(a), tanh(a)."""

    expr = expand_hyperbolic_trig(expr)
    return approximate_exp(expr)

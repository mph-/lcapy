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


def _approximate_exp_pade(expr, order=1, numer_order=None):

    if numer_order is None:
        numer_order = order

    def query(expr):
        return expr.is_Function and expr.func == exp

    def value(expr):
        arg = expr.args[0]

        if order == 1 and numer_order == 1:
            return (2 + arg) / (2 - arg)
        elif order == 2 and numer_order == 2:
            return (12 + 6 * arg + arg**2) / (12 - 6 * arg + arg**2)

        from math import factorial

        m = numer_order
        n = order
        numer = 0
        denom = 0

        for k in range(m + 1):
            scale = (factorial(m + n - k) * factorial(m)
                     ) // (factorial(k) * factorial(m - k) * factorial(m))
            numer += scale * arg**k

        for k in range(n + 1):
            scale = (factorial(m + n - k) * factorial(n)
                     ) // (factorial(k) * factorial(n - k) * factorial(m))
            denom += scale * (-arg)**k

        return numer / denom

    return expr.replace(query, value)


def approximate_exp(expr, method='pade', order=1, numer_order=None):
    """Approximate exp(a).  The best time-domain response (without a jump)
    is achieved with 'numer_order == order - 1'.  The best
    frequency-domain response is achieved with numer_order ==
    order."""

    if method != 'pade':
        raise ValueError('Method %s unsupported, must be pade')

    return _approximate_exp_pade(expr, order, numer_order)


def approximate_hyperbolic_trig(expr, method='pade', order=1, numer_order=None):
    """Approximate cosh(a), sinh(a), tanh(a)."""

    expr = expand_hyperbolic_trig(expr)
    return approximate_exp(expr, method=method, order=order,
                           numer_order=numer_order)


def approximate_dominant_terms(expr, defs, threshold=0.01):

    terms = expr.as_ordered_terms()

    # Cannot substitute first since SymPy reorders the args.

    absvalmax = None
    for term in terms:
        nterm = term.subs(defs)
        if not nterm.is_constant():
            continue
        absval = abs(nterm)
        if absvalmax is None or absval > absvalmax:
            absvalmax = absval

    newexpr = 0
    for term in terms:
        nterm = term.subs(defs)
        absval = abs(nterm)
        if (not nterm.is_constant()) or (absval > threshold * absvalmax):
            newexpr += term

    return newexpr


def approximate_dominant(expr, defs, threshold=0.01):
    """Approximate expression using ball-park numerical values for the
    symbols to decide which terms in a sum dominate the sum.  A term
    is neglected if its absolute value is below `threshold` times the
    maximum absolute value of the terms in the sum.

    `defs` is a dict, set, list, or tuple of symbol definitions."""

    def query(expr):

        return expr.is_Add

    def value(expr):
        return approximate_dominant_terms(expr, defs, threshold)

    return expr.replace(query, value)


def approximate_order(expr, var, order):
    """Approximate expression by reducing order of polynomial to
    specified order."""

    from sympy import O

    return (expr.expand() + O(var ** (order + 1))).removeO()

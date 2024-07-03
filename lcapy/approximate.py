"""This module contains functions for approximating SymPy expressions.

Copyright 2021--2024 Michael Hayes, UCECE

"""

from .simplify import expand_hyperbolic_trig
from sympy import exp, diff, eye, zeros, Matrix, Poly
from math import factorial


def approximate_fractional_power(expr, var, ndegree=2, ddegree=2):
    """This is an experimental method to approximate s**a, where a is
    fractional, with a rational function using a Pade approximant
    of order [ndegree, ddegree]."""

    v = var

    def query(expr):

        if not expr.is_Pow:
            return False
        # TODO handle (b * v)**a
        if expr.args[0] != v:
            return False
        if expr.args[1].is_Number and not expr.args[1].is_Integer:
            return True
        if expr.args[1].is_Symbol and not expr.args[1].is_Integer:
            return True
        if (expr.args[1].is_Mul and expr.args[1].args[0] == -1
            and expr.args[1].args[1].is_Symbol
            and not expr.args[1].args[1].is_Integer):
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

    if ndegree != ddegree:
        raise ValueError('Require ndegree == ddegree')

    if ndegree == 1:
        value = value1
    elif ndegree == 2:
        value = value2
    else:
        raise ValueError('Can only handle degree 1 and 2 at the moment')

    expr = expr.replace(query, value)

    return expr


def approximate_exp(expr, ndegree=1, ddegree=1):
    """Approximate exp(a) with Pade approximant.  The best time-domain
    response (without a jump) is achieved with 'ndegree == ddegree -
    1'.  The best frequency-domain response is achieved with ndegree
    == ddegree.

    """

    def query(expr):
        return expr.is_Function and expr.func == exp

    def value(expr):
        arg = expr.args[0]

        if ddegree == 1 and ndegree == 1:
            return (2 + arg) / (2 - arg)
        elif ddegree == 2 and ndegree == 2:
            return (12 + 6 * arg + arg**2) / (12 - 6 * arg + arg**2)

        from math import factorial

        m = ndegree
        n = ddegree
        numer = 0
        denom = 0

        for k in range(m + 1):
            scale = (factorial(m + n - k) * factorial(m)) \
                // (factorial(k) * factorial(m - k) * factorial(m))
            numer += scale * arg**k

        for k in range(n + 1):
            scale = (factorial(m + n - k) * factorial(n)) \
                // (factorial(k) * factorial(n - k) * factorial(m))
            denom += scale * (-arg)**k

        return numer / denom

    return expr.replace(query, value)


def approximate_hyperbolic_trig(expr, ndegree=1, ddegree=1):
    """Approximate cosh(a), sinh(a), tanh(a)."""

    expr = expand_hyperbolic_trig(expr)
    return approximate_exp(expr, ndegree=ndegree, ddegree=ddegree)


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


def approximate_degree(expr, var, degree):
    """Approximate expression by reducing degree of polynomial to
    specified degree."""

    from sympy import O

    return (expr.expand() + O(var ** (degree + 1))).removeO()


def approximate_taylor_coeffs(expr, var, degree=1, var0=0):

    coeffs = [expr.subs(var, var0)]

    prev = expr

    for n in range(1, degree + 1):

        prev = diff(prev, var)

        coeffs.append(prev.subs(var, var0) / factorial(n))

    return coeffs


def approximate_taylor(expr, var, degree=1, var0=0):
    """Approximate expression using a Taylor series
    around `var = var0` to degree `degree`."""

    result = expr.subs(var, var0)

    prev = expr

    for n in range(1, degree + 1):

        prev = diff(prev, var)

        result += prev.subs(var, var0) * \
            (var - var0)**n / factorial(n)

    return result


def approximate_pade_coeffs_from_taylor_coeffs(a, m):

    N = len(a) - 1
    n = N - m
    if n < 0:
        raise ValueError('Degree %d must be smaller than %d' % (m, N))
    A = eye(N + 1, n + 1)
    B = zeros(N + 1, m)

    for col in range(m):
        for row in range(col + 1, N + 1):
            B[row, col] = -a[row - col - 1]

    C = A.hstack(A, B)
    pq = C.solve(Matrix(a))
    p = pq[:n + 1]

    q = [1]
    for col in range(1, m + 1):
        q.append(pq[n + col])

    return p[::-1], q[::-1]


def approximate_pade_coeffs(expr, var, ndegree=2, ddegree=2, var0=0):
    """Approximate expression using a Pade series with numerator degree
    `ndegree` and denominator degree `ddegree`."""

    coeffs = approximate_taylor_coeffs(expr, var, ndegree + ddegree,
                                       var0=var0)

    return approximate_pade_coeffs_from_taylor_coeffs(coeffs, ddegree)


def approximate_pade(expr, var, ndegree=2, ddegree=2, var0=0):

    p, q = approximate_pade_coeffs(expr, var, ndegree, ddegree, var0)

    return Poly.from_list(p, gens=var - var0) \
        / Poly.from_list(q, gens=var - var0)

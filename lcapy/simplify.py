"""This module contains functions for simplifying expressions.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from sympy import Add, Mul, DiracDelta, Heaviside, Integral, re, im
from sympy import oo, sin, cos, sqrt, atan2, pi, Symbol, solve, Min, Max
from sympy import cosh, sinh, tanh, exp
from .extrafunctions import UnitStep, UnitImpulse, rect, dtrect


def simplify_dirac_delta_product_term(expr):
    """ Simplify f(t) * DiracDelta(t) to f(0) * DiracDelta(t)."""

    if not expr.has(DiracDelta):
        return expr

    def query(expr):

        return expr.is_Mul and expr.has(DiracDelta)

    def value(expr):

        arg = None
        factors = expr.args
        dirac = None
        parts = []
        for factor in factors:
            if factor.is_Function and factor.func == DiracDelta:
                arg = factor.args[0]
                dirac = factor
            else:
                parts.append(factor)

        if arg is None or not arg.has(Symbol):
            return expr
        results = solve(arg, dict=True)
        # Note, the eval method is called for functions.
        const = Mul(*parts).subs(results[0])
        return const * dirac

    return expr.replace(query, value)


def simplify_dirac_delta_product(expr, expand=False):
    """Simplify f(t) * DiracDelta(t) to f(0) * DiracDelta(t)."""

    if not expr.has(DiracDelta):
        return expr

    # Could also convert delta(a * t) to delta(t) / a

    if not expand:
        return simplify_dirac_delta_product_term(expr)

    terms = expr.expand().as_ordered_terms()

    return Add(*[simplify_dirac_delta_product_term(term) for term in terms])


def simplify_dirac_delta(expr, var=None):

    if not expr.has(DiracDelta):
        return expr

    expr = simplify_dirac_delta_product(expr)
    if var is not None:

        # Convert delta(a * t) to delta(t) / a
        expr = expr.expand(diracdelta=True, wrt=var)
    return expr


# def simplify_heaviside_product(expr):
#
#     heaviside_products = []
#
#     def pre(expr):
#         if (expr.is_Mul and expr.args[0].func == Heaviside and
#               expr.args[1].func == Heaviside):
#             heaviside_products.append(expr)
#
#         for arg in expr.args:
#             pre(arg)
#
#     pre(expr)
#
#     for product in heaviside_products:
#         # TODO
#         pass
#
#     return expr


def simplify_power(expr):

    powers = []

    def pre(expr):
        if (expr.is_Pow and expr.args[0].func in (Heaviside, UnitStep, rect, dtrect) and
                expr.args[1].is_constant):
            powers.append(expr)

        for arg in expr.args:
            pre(arg)

    pre(expr)

    for power in powers:
        expr = expr.replace(power, power.args[0])
    return expr


def simplify_heaviside_integral(expr):

    if not expr.has(Integral):
        return expr

    def query(expr):

        if not isinstance(expr, Integral):
            return False
        return expr.has(Heaviside) or expr.has(UnitStep)

    def value(expr):

        integrand = expr.args[0]
        var = expr.args[1][0]
        lower_limit = expr.args[1][1]
        upper_limit = expr.args[1][2]

        # Rewrite integral limits if Heaviside is a factor of the
        # integrand.

        result = 1
        for factor in integrand.as_ordered_factors():
            if isinstance(factor, (Heaviside, UnitStep)):
                arg = factor.args[0]
                if arg == var:
                    lower_limit = Max(lower_limit, 0)
                    factor = 1
                elif arg == -var:
                    upper_limit = Min(upper_limit, 0)
                    factor = 1
                elif (arg.is_Add and arg.args[1].is_Mul and
                      arg.args[1].args[0] == -1 and arg.args[1].args[1] == var):
                    upper_limit = Min(upper_limit, arg.args[0])
                    # Cannot remove Heaviside function in general.

            result *= factor

        ret = Integral(result, (var, lower_limit, upper_limit))
        return ret

    expr = expr.replace(query, value)

    return expr


def simplify_heaviside_scale(expr, var):

    terms = expr.as_ordered_terms()
    if len(terms) > 1:
        result = 0
        for term in terms:
            result += simplify_heaviside_scale(term, var)
        return result

    def query(expr):

        return expr.is_Function and expr.func in (Heaviside, UnitStep)

    def value(expr):

        arg = expr.args[0]

        if not arg.as_poly(var).is_linear:
            return expr

        arg = arg.expand()
        a = arg.coeff(var, 1)
        b = arg.coeff(var, 0)
        if a == 0:
            return expr

        return expr.func(var + (b / a).cancel())

    return expr.replace(query, value)


def simplify_heaviside(expr, var=None):

    if not expr.has(Heaviside) and not expr.has(UnitStep):
        return expr

    expr = simplify_heaviside_integral(expr)
    expr = simplify_power(expr)
    if var is not None:
        expr = simplify_heaviside_scale(expr, var)
    return expr


def simplify_rect(expr, var=None):

    if not expr.has(rect) and not expr.has(dtrect):
        return expr

    expr = simplify_power(expr)
    return expr


def simplify_sin_cos(expr, as_cos=False, as_sin=False):

    if not (expr.has(sin) and expr.has(cos)):
        return expr

    terms = expr.expand().as_ordered_terms()

    rest = 0
    cos_parts = []
    sin_parts = []

    for term in terms:
        if term.has(sin):
            sin_parts.append(term)
        elif term.has(cos):
            cos_parts.append(term)
        else:
            rest += term

    if cos_parts == [] or sin_parts == []:
        return expr

    result = 0

    # Use list to make copy for iterator.
    for cos_part in list(cos_parts):
        cfactors = cos_part.expand().as_ordered_factors()

        for sin_part in list(sin_parts):
            sfactors = sin_part.expand().as_ordered_factors()

            commonfactors = []
            for factor in cfactors:
                if factor in sfactors:
                    commonfactors.append(factor)

            for factor in commonfactors:
                sfactors.remove(factor)
                cfactors.remove(factor)

            cos_factor = None
            sin_factor = None
            for cfactor in cfactors:
                if cfactor.has(cos):
                    cos_factor = cfactor
                    break

            for sfactor in sfactors:
                if sfactor.has(sin):
                    sin_factor = sfactor
                    break

            if cos_factor is None or sin_factor is None:
                continue

            if cos_factor.args[0] != sin_factor.args[0]:
                continue

            cfactors.remove(cos_factor)
            sfactors.remove(sin_factor)

            c = Mul(*cfactors)
            s = Mul(*sfactors)
            A = sqrt(c * c + s * s) * Mul(*commonfactors)
            phi = atan2(s, c)

            if as_sin:
                result += A * \
                    sin(cos_factor.args[0] - phi + pi / 2, evaluate=False)

            if as_cos:
                result += A * cos(cos_factor.args[0] - phi, evaluate=False)

            # SymPy will choose sin or cos as convenient.
            if not as_sin and not as_cos:
                result += A * cos(cos_factor.args[0] - phi)

            cos_parts.remove(cos_part)
            sin_parts.remove(sin_part)
            break

    result = Add(result, rest, *cos_parts, *sin_parts)
    return result


def simplify_unit_impulse(expr, var=None):

    if not expr.has(UnitImpulse):
        return expr

    def query(expr):

        return expr.is_Function and expr.func is UnitImpulse

    def value(expr):

        arg = expr.args[0]

        if not arg.as_poly(var).is_linear:
            return expr

        arg = arg.expand()
        a = arg.coeff(var, 1)
        b = arg.coeff(var, 0)
        if a == 0:
            return expr

        return expr.func(var + (b / a).cancel())

    return expr.replace(query, value)


def simplify_conjugates(expr):

    factors = expr.as_ordered_factors()

    if len(factors) > 1:
        sfactors = []
        for factor in factors:
            sfactors.append(simplify_conjugates(factor))
        return Mul(*sfactors)

    terms = expr.expand().as_ordered_terms()

    sterms = []
    for m, term in enumerate(terms):
        if term == 0:
            continue
        cterm = term.conjugate()

        for m1, term1 in enumerate(terms[m + 1:]):
            if cterm == term1:
                terms[m1 + m + 1] = 0
                term = 2 * re(term)
                break
            elif cterm == -term1:
                terms[m1 + m + 1] = 0
                term = 2 * im(term)
                break

        sterms.append(term)

    return Add(*sterms)


def expand_hyperbolic_trig(expr):

    # Note, rewrite(exp) does this except it also converts s**2 to
    # exp(2 * log(s)).

    def query(expr):
        return expr.is_Function and expr.func in (cosh, sinh, tanh)

    def value(expr):
        arg = expr.args[0]

        if expr.func == cosh:
            return exp(arg) + exp(-arg)
        elif expr.func == sinh:
            return exp(arg) - exp(-arg)
        elif expr.func == tanh:
            return (exp(arg) - exp(-arg)) / (exp(arg) + exp(-arg))
        else:
            raise RuntimeError('Internal error')

    return expr.replace(query, value)


# Could have simplify_sum for things like Sum(0**m * x(m), (m, 0, oo))
# that simplifies to x(0).

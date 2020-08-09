"""This module contains functions for simplifying expressions.

Copyright 2020 Michael Hayes, UCECE

"""

from sympy import Add, Mul, DiracDelta, Heaviside, Integral, oo


def simplify_dirac_delta_term(expr):
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

        # TODO: handle DiracDelta(symbol + constant).
        if arg is None or not arg.is_Symbol:
            return expr
        const = Mul(*parts).subs(arg, 0)
        return const * dirac

    return expr.replace(query, value)


def simplify_dirac_delta(expr, expand=False):
    """Simplify f(t) * DiracDelta(t) to f(0) * DiracDelta(t)."""

    if not expr.has(DiracDelta):
        return expr

    if not expand:
        return simplify_dirac_delta_term(expr)
    
    terms = expr.expand().as_ordered_terms()

    return Add(*[simplify_dirac_delta_term(term) for term in terms])


def simplify_heaviside(expr):

    if not expr.has(Heaviside):
        return expr

    # H(t) x H(t) = H(t)

    if not expr.has(Integral):
        return expr

    def query(expr):

        if not isinstance(expr, Integral):
            return False
        return expr.has(Heaviside)

    def value(expr):

        integrand = expr.args[0]
        var = expr.args[1][0]
        lower_limit = expr.args[1][1]
        upper_limit = expr.args[1][2]

        # Rewrite integral limits if Heaviside is a factor of the
        # integrand.
        # TODO, be more clever if limits not infinite, say using min(l1, l2)
        
        result = 1
        for factor in integrand.as_ordered_factors():
            if isinstance(factor, Heaviside):
                arg = factor.args[0]
                if arg == var and lower_limit == -oo:
                    lower_limit = 0
                    factor = 1
                elif arg == -var and upper_limit == oo:
                    upper_limit = 0
                    factor = 1
                elif arg.is_Add and arg.has(-var) and upper_limit == -oo:
                    upper_limit = var + arg
                    factor = 1
                elif arg.is_Add and arg.has(var) and lower_limit == oo:
                    lower_limit = var - arg
                    factor = 1                                        
                    
            result *= factor

        return Integral(result, (var, lower_limit, upper_limit))
    
    expr = expr.replace(query, value)
    
    return expr

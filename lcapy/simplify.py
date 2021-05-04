"""This module contains functions for simplifying expressions.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from sympy import Add, Mul, DiracDelta, Heaviside, Integral, oo, sin, cos, sqrt, atan2, pi, Symbol, solve


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
        if arg is None or not arg.has(Symbol):
            return expr
        results = solve(arg, dict=True)
        # Note, the eval method is called for functions.
        const = Mul(*parts).subs(results[0])
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


def simplify_heaviside_power(expr):

    heaviside_powers = []    
    
    def pre(expr):
        if (expr.is_Pow and expr.args[0].func == Heaviside and
              expr.args[1].is_constant):
            heaviside_powers.append(expr)            
        
        for arg in expr.args:
            pre(arg)

    pre(expr)

    for power in heaviside_powers:
        expr = expr.replace(power, power.args[0])            
    return expr


def simplify_heaviside_integral(expr):

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


def simplify_heaviside(expr):

    if not expr.has(Heaviside):
        return expr
    
    expr = simplify_heaviside_integral(expr)
    expr = simplify_heaviside_power(expr)    
    return expr
    

def simplify_sin_cos(expr, as_cos=False, as_sin=False):

    if not (expr.has(sin) and expr.has(cos)):
        return expr
    
    terms = expr.expand().as_ordered_terms()

    rest = 0
    cos_part = None
    sin_part = None    
    
    for term in terms:
        if term.has(sin) and sin_part is None:
            sin_part = term
        elif term.has(cos) and cos_part is None:
            cos_part = term
        else:
            rest += term

    if cos_part is None or sin_part is None:
        return expr

    cfactors = cos_part.expand().as_ordered_factors()
    sfactors = sin_part.expand().as_ordered_factors()

    commonfactors = []
    for factor in cfactors:
        if factor in sfactors:
            commonfactors.append(factor)

    for factor in commonfactors:
        sfactors.remove(factor)
        cfactors.remove(factor)

    cosfactor = None
    sinfactor = None    
    for cfactor in cfactors:
        if cfactor.has(cos):
            cosfactor = cfactor
            break
        
    for sfactor in sfactors:
        if sfactor.has(sin):
            sinfactor = sfactor
            break
        
    if cosfactor is None or sinfactor is None:
        return expr

    if cosfactor.args[0] != sinfactor.args[0]:
        return expr
        
    cfactors.remove(cosfactor)
    sfactors.remove(sinfactor)    

    c = Mul(*cfactors)
    s = Mul(*sfactors)
    A = sqrt(c * c + s * s) * Mul(*commonfactors)
    phi = atan2(s, c)

    if as_sin:
        return rest + A * sin(cosfactor.args[0] - phi + pi / 2, evaluate=False)

    if as_cos:
        return rest + A * cos(cosfactor.args[0] - phi, evaluate=False)

    # SymPy will choose sin or cos as convenient.
    return rest + A * cos(cosfactor.args[0] - phi)

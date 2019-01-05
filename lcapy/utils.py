import sympy as sym

def factor_const(expr, t):

    # Perhaps use expr.as_coeff_Mul() ?

    rest = sym.S.One
    const = sym.S.One
    for factor in expr.as_ordered_factors():
        # Cannot use factor.is_constant() since SymPy 1.2, 1.3
        # barfs for Heaviside(t) and DiracDelta(t)
        if not factor.has(t):
            const *= factor
        else:
            rest *= factor
    return const, rest


def scale_shift(expr, t):

    if not expr.has(t):
        raise ValueError('Expression does not contain %s: %s' % (t, expr))

    terms = expr.as_ordered_terms()
    if len(terms) > 2:
        raise ValueError('Expression has too many terms: %s' % expr)

    if len(terms) == 1:
        return terms[0] / t, sym.S.Zero

    scale = terms[0] / t
    if not scale.is_constant():
        raise ValueError('Expression not a scale and shift: %s' % expr)

    return scale, terms[1]

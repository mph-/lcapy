"""This module provides support for Fourier transforms.  It acts as a
wrapper for SymPy's Fourier transform.  It calculates the bilateral
Fourier transform using:

   S(f) = \int_{-\infty}^{\infty} s(t) e^{-j * 2* \pi * t} dt

It also allows functions that strictly do not have a Fourier transform
by using Dirac deltas.  For example, a, cos(a * t), sin(a * t), exp(j
* a * t).


Copyright 2016 Michael Hayes, UCECE

"""

import sympy as sym

def fourier_sympy(expr, t, f):

    result = sym.fourier_transform(expr, t, f)
    if expr != 0 and result == 0:
        # There is a bug in SymPy where it returns 0.
        raise ValueError('Could not compute Fourier transform for ' + str(expr))

    return result


def fourier_term(expr, t, f):

    var = sym.Symbol(str(t))
    expr = expr.replace(var, t)

    if expr.has(sym.function.AppliedUndef) and expr.args[0] == t:
        # TODO, handle things like 3 * v(t), a * v(t), 3 * t * v(t), v(t-T),
        # v(4 * a * t), etc.
        if not isinstance(expr, sym.function.AppliedUndef):
            raise ValueError('Could not compute Fourier transform for ' + str(expr))

        # Convert v(t) to V(f), etc.
        name = expr.func.__name__
        name = name[0].upper() + name[1:] + '(f)'
        return sym.sympify(name)

    # Check for constant.
    if not expr.has(t):
        return expr * sym.DiracDelta(f)

    one = sym.sympify(1)
    const = one
    other = one
    exps = one
    factors = expr.as_ordered_factors()    
    for factor in factors:
        if not factor.has(t):
            const *= factor
        else:
            if factor.is_Function and factor.func == sym.exp:
                exps *= factor
            else:
                other *= factor

    if other != 1:
        if other == t:
            return const * sym.I * 2 * sym.pi * sym.DiracDelta(f, 1)
        if other == t**2:
            return const * (sym.I * 2 * sym.pi)**2 * sym.DiracDelta(f, 2)

        # Punt and use SymPy.  Should check for t**n, t**n * exp(-a * t), etc.
        return fourier_sympy(expr, t, f)

    args = exps.args[0]
    foo = args / t
    if foo.has(t):
        # Have exp(a * t**n), SymPy might be able to handle this
        return fourier_sympy(expr, t, f)

    if exps != 1:
        return const * sym.DiracDelta(f - foo / (-sym.I * 2 * sym.pi))
        
    return fourier_sympy(expr, t, f)


def fourier_transform(expr, t, f):
    """Compute bilateral Fourier transform of expr.

    Undefined functions such as v(t) are converted to V(f)

    This also handles some expressions that do not really have a Fourier
    transform, such as a, cos(a*t), sin(a*t), exp(I * a * t).

    """

    # The variable may have been created with different attributes,
    # say when using sym.sympify('DiracDelta(t)') since this will
    # default to assuming that t is complex.  So if the symbol has the
    # same representation, convert to the desired one.

    var = sym.Symbol(str(t))
    expr = expr.replace(var, t)

    orig_expr = expr

    if expr.has(sym.cos) or expr.has(sym.sin):
        expr = expr.rewrite(sym.exp)

    terms = expr.expand().as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += fourier_term(term, t, f)
    except ValueError:
        raise ValueError('Could not compute Fourier transform for ' + str(orig_expr))

    return result


def test():

     t, f, a = sym.symbols('t f a', real=True)

     print(fourier_transform(a, t, f))
     print(fourier_transform(sym.exp(-sym.I * 2 * sym.pi * a * t), t, f))
     print(fourier_transform(sym.cos(2 * sym.pi * a * t), t, f))
     print(fourier_transform(sym.sin(2 * sym.pi * a * t), t, f))
     print(fourier_transform(a * t, t, f))
     print(fourier_transform(sym.exp(-a * t) * sym.Heaviside(t), t, f))



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

fourier_cache = {}

def fourier_sympy(expr, t, f):

    result = sym.fourier_transform(expr, t, f)
    if expr != 0 and result == 0:
        # There is a bug in SymPy where it returns 0.
        raise ValueError('Could not compute Fourier transform for ' + str(expr))

    return result


def fourier_term(expr, t, f, inverse=False):

    if expr.has(sym.function.AppliedUndef) and expr.args[0] == t:
        # TODO, handle things like 3 * v(t), a * v(t), 3 * t * v(t), v(t-T),
        # v(4 * a * t), etc.
        if not isinstance(expr, sym.function.AppliedUndef):
            raise ValueError('Could not compute Fourier transform for ' + str(expr))

        # Convert v(t) to V(f), etc.
        name = expr.func.__name__
        if inverse:
            name = name[0].lower() + name[1:] + '(%s)' % -f
        else:
            name = name[0].upper() + name[1:] + '(%s)' % f
        return sym.sympify(name)

    # Check for constant.
    if not expr.has(t):
        if f.is_Mul and f.args[0] < 0:
            f = -f
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

    if other != 1 and exps == 1:
        if other == t:
            return const * sym.I * 2 * sym.pi * sym.DiracDelta(f, 1)
        if other == t**2:
            return const * (sym.I * 2 * sym.pi)**2 * sym.DiracDelta(f, 2)

        # Sympy incorrectly gives exp(-a * t) instead of exp(-a * t) *
        # Heaviside(t)
        if other.is_Pow and other.args[1] == -1:
            foo = other.args[0]
            if foo.is_Add and foo.args[1].has(t):
                bar = foo.args[1] / t
                if not bar.has(t) and bar.has(sym.I):
                    a = -(foo.args[0] * 2 * sym.pi * sym.I) / bar
                    return const * sym.exp(-a * f) * sym.Heaviside(f * sym.sign(a))

        # Punt and use SymPy.  Should check for t**n, t**n * exp(-a * t), etc.
        return fourier_sympy(expr, t, f)

    args = exps.args[0]
    foo = args / t
    if foo.has(t):
        # Have exp(a * t**n), SymPy might be able to handle this
        return fourier_sympy(expr, t, f)

    if exps != 1 and foo.has(sym.I):
        return const * sym.DiracDelta(f - foo / (sym.I * 2 * sym.pi))
        
    return fourier_sympy(expr, t, f)


def fourier_transform(expr, t, f, inverse=False):
    """Compute bilateral Fourier transform of expr.

    Undefined functions such as v(t) are converted to V(f)

    This also handles some expressions that do not really have a Fourier
    transform, such as a, cos(a * t), sin(a * t), exp(I * a * t).

    """

    key = (expr, t, f)
    if key in fourier_cache:
        return fourier_cache[key]
    
    # Hack for debugging.  Otherwise sym.sympify will convert Expr
    # types to string and then re-parse.  Unfortunately, we change I
    # to j when printing and so j gets converted into a symbol and not
    # the imaginary unit.
    if hasattr(expr, 'expr'):
        expr = expr.expr
    if hasattr(t, 'expr'):
        t = t.expr
    if hasattr(f, 'expr'):
        f = f.expr

    expr = sym.sympify(expr)
    t = sym.sympify(t)
    f = sym.sympify(f)

    if inverse:
        t, f = f, -t

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
            result += fourier_term(term, t, f, inverse=inverse)
    except ValueError:
        raise ValueError('Could not compute Fourier transform for ' + str(orig_expr))

    fourier_cache[key] = result
    return result



def inverse_fourier_transform(expr, f, t):
    """Compute bilateral inverse Fourier transform of expr.

    Undefined functions such as V(f) are converted to v(t)

    This also handles some expressions that do not really have an
    inverse Fourier transform, such as a, cos(a * f), sin(a * f), exp(I *
    a * f).

    """

    result = fourier_transform(expr, t, f, inverse=True)
    return sym.simplify(result)


def test():

     t, f, a = sym.symbols('t f a', real=True)
     a = sym.symbols('a', positive=True)

     print(fourier_transform(a, t, f))
     print(fourier_transform(sym.exp(-sym.I * 2 * sym.pi * a * t), t, f))
     print(fourier_transform(sym.cos(2 * sym.pi * a * t), t, f))
     print(fourier_transform(sym.sin(2 * sym.pi * a * t), t, f))
     print(fourier_transform(a * t, t, f))
     print(fourier_transform(sym.exp(-a * t) * sym.Heaviside(t), t, f))
     print(inverse_fourier_transform(a, f, t))
     print(inverse_fourier_transform(1 / (sym.I * 2 * sym.pi * f + a), f, t))

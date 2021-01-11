"""This module provides support for Fourier transforms.  It acts as a
quantity for SymPy's Fourier transform.  It calculates the bilateral
Fourier transform using:

   S(f) = \int_{-\infty}^{\infty} s(t) e^{-j * 2 * \pi * t} dt

It also allows functions that strictly do not have a Fourier transform
by using Dirac deltas.  For example, a, cos(a * t), sin(a * t), exp(j
* a * t).


Copyright 2016--2021 Michael Hayes, UCECE

"""

# TODO:
# Add DiracDelta(t, n)
# Simplify  (-j * DiracDelta(f - 1) +j * DiracDelta(f + 1)).inverse_fourier()
# This should give 2 * sin(2 * pi * t)

import sympy as sym
from .sym import sympify, AppliedUndef, j, pi, symsimplify
from .utils import factor_const, scale_shift
from .extrafunctions import rect, sinc

__all__ = ('FT', 'IFT')


fourier_cache = {}

def fourier_sympy(expr, t, f):

    result = sym.fourier_transform(expr, t, f)
    if expr != 0 and result == 0:
        # There is a bug in SymPy where it returns 0.
        raise ValueError('Could not compute Fourier transform for ' + str(expr))

    if isinstance(result, sym.FourierTransform):
        raise ValueError('Could not compute Fourier transform for ' + str(expr))
    
    return result


def fourier_func(expr, t, f, inverse=False):

    if not isinstance(expr, AppliedUndef):
        raise ValueError('Expecting function for %s' % expr)

    scale, shift = scale_shift(expr.args[0], t)

    fsym = sympify(str(f))
    
    # Convert v(t) to V(f), etc.
    name = expr.func.__name__
    if inverse:
        func = name[0].lower() + name[1:] + '(%s)' % f
    else:
        func = name[0].upper() + name[1:] + '(%s)' % f

    result = sympify(func).subs(fsym, f / scale) / abs(scale)

    if shift != 0:
        if inverse:
            shift = -shift
        result = result * sym.exp(2 * sym.I * sym.pi * f * shift / scale)
    
    return result


def fourier_function(expr, t, f, inverse=False):

    # Handle expressions with a function of FOO, e.g.,
    # v(t), v(t) * y(t),  3 * v(t) / t, v(4 * a * t), etc.,
    
    if not expr.has(AppliedUndef):
        raise ValueError('Could not compute Fourier transform for ' + str(expr))

    const, expr = factor_const(expr, t)
    
    if isinstance(expr, AppliedUndef):
        return fourier_func(expr, t, f, inverse) * const
    
    tsym = sympify(str(t))
    expr = expr.subs(tsym, t)

    rest = sym.S.One
    undefs = []
    for factor in expr.as_ordered_factors():
        if isinstance(factor, AppliedUndef):
            if factor.args[0] != t:
                raise ValueError('Weird function %s not of %s' % (factor, t))
            undefs.append(factor)
        else:
            rest *= factor

    if rest.has(AppliedUndef):
        # Have something like 1/v(t)
        raise ValueError('Cannot compute Fourier transform of %s' % rest)
            
    exprs = undefs
    if rest.has(t):
        exprs = exprs + [rest]
        rest = sym.S.One

    result = fourier_term(exprs[0], t, f, inverse) * rest
        
    if len(exprs) == 1:
        return result * const

    dummy = 'tau' if inverse else 'nu'

    for m in range(len(exprs) - 1):
        if m == 0:
            nu = sympify(dummy)
        else:
            nu = sympify(dummy + '_%d' % m)
        expr2 = fourier_term(exprs[m + 1], t, f, inverse)
        result = sym.Integral(result.subs(f, f - nu) * expr2.subs(f, nu),
                              (nu, -sym.oo, sym.oo))
    
    return result * const

def fourier_term(expr, t, f, inverse=False):

    const, expr = factor_const(expr, t)
    
    if isinstance(expr, AppliedUndef):
        return fourier_func(expr, t, f, inverse) * const
    
    # TODO add u(t) <-->  delta(f) / 2 - j / (2 * pi * f)
    
    if expr.has(AppliedUndef):
        # Handle v(t), v(t) * y(t),  3 * v(t) / t etc.
        return fourier_function(expr, t, f, inverse) * const

    # Check for constant.
    if not expr.has(t):
        return expr * sym.DiracDelta(f) * const

    one = sym.S.One
    const1 = const
    other = one
    exps = one
    factors = expr.expand().as_ordered_factors()    
    for factor in factors:
        if not factor.has(t):
            const1 *= factor
        else:
            if factor.is_Function and factor.func == sym.exp:
                exps *= factor
            else:
                other *= factor

    sf = -f if inverse else f
    
    if other != 1 and exps == 1:
        if other == t:
            return const1 * sym.I / (2 * sym.pi) * sym.DiracDelta(sf, 1)
        elif other == t**2:
            return -const1 / (2 * sym.pi)**2 * sym.DiracDelta(sf, 2)
        # TODO check for other powers of t...
        elif other == sym.sign(t):
            return const1 / (sym.I * sym.pi * sf)
        elif other == sym.sign(t) * t:
            return -const1 * 2 / (2 * sym.pi * f)**2
        elif other == sym.Heaviside(t):
            return const1 / (sym.I * 2 * sym.pi * f) + const1 * sym.DiracDelta(sf) / 2
        elif other == 1 / t:
            return -const1 * sym.I * sym.pi * sym.sign(sf)
        elif other == 1 / t**2:
            return -const1 * 2 * sym.pi**2 * sf * sym.sign(sf)        
        elif other.is_Function and other.func == sym.Heaviside and other.args[0].has(t):
            # TODO, generalise use of similarity and shift theorems for other functions and expressions
            scale, shift = scale_shift(other.args[0], t)
            return (const1 / (sym.I * 2 * sym.pi * sf / scale) / abs(scale) + const1 * sym.DiracDelta(sf) / 2) * sym.exp(sym.I * 2 * sym.pi * sf /scale * shift)
        elif other == sym.Heaviside(t) * t:
            return -const1 / (2 * sym.pi * f)**2 + const1 * sym.I * sym.DiracDelta(sf, 1) / (4 * pi)
        elif other.is_Function and other.func == sinc and other.args[0].has(t):
            scale, shift = scale_shift(other.args[0], t)            
            return const1 * rect(sf / scale) * sym.exp(sym.I * 2 * sym.pi * sf /scale * shift) / abs(scale)
        elif other.is_Function and other.func == rect and other.args[0].has(t):        
            return const1 * sinc(sf / scale) * sym.exp(sym.I * 2 * sym.pi * sf /scale * shift) / abs(scale)

        # Sympy incorrectly gives exp(-a * t) instead of exp(-a * t) *
        # Heaviside(t)
        if other.is_Pow and other.args[1] == -1:
            foo = other.args[0]
            if foo.is_Add and foo.args[1].has(t):
                bar = foo.args[1] / t
                if not bar.has(t) and bar.has(sym.I):
                    a = -(foo.args[0] * 2 * sym.pi * sym.I) / bar
                    return const1 * sym.exp(-a * sf) * sym.Heaviside(sf * sym.sign(a))

        if expr == t * sym.DiracDelta(t, 1):
            return const * sf / (-sym.I * 2 * sym.pi)
                
        # Punt and use SymPy.  Should check for t**n, t**n * exp(-a * t), etc.
        return const * fourier_sympy(expr, t, sf)

    args = exps.args[0]
    foo = args / t
    if foo.has(t):
        # Have exp(a * t**n), SymPy might be able to handle this
        return const * fourier_sympy(expr, t, sf)

    if exps != 1 and foo.has(sym.I):
        return const1 * sym.DiracDelta(sf - foo / (sym.I * 2 * sym.pi))
        
    return const * fourier_sympy(expr, t, sf)


def fourier_transform(expr, t, f, inverse=False, evaluate=True):
    """Compute bilateral Fourier transform of expr.

    Undefined functions such as v(t) are converted to V(f)

    This also handles some expressions that do not really have a
    Fourier transform, such as a, cos(a * t), sin(a * t), exp(I * a *
    t).  These expressions all require the use of the Dirac delta."""

    if expr.is_Equality:
        return sym.Eq(fourier_transform(expr.args[0], t, f, inverse),
                      fourier_transform(expr.args[1], t, f, inverse))

    if not evaluate:
        if inverse:
            result = sym.Integral(expr * sym.exp(j * 2 * pi * f * t), (f, -sym.oo, sym.oo))
        else:
            result = sym.Integral(expr * sym.exp(-j * 2 * pi * f * t), (t, -sym.oo, sym.oo))
        return result

    key = (expr, t, f, inverse)
    if key in fourier_cache:
        return fourier_cache[key]

    if not inverse and expr.has(f):
        raise ValueError('Cannot Fourier transform for expression %s that depends on %s' % (expr, f))


    if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
        raise ValueError('Cannot Fourier transform for expression %s that is unknown for t < 0' % expr)
    
    if inverse:
        t, f = f, t

    # The variable may have been created with different attributes,
    # say when using sym.sympify('DiracDelta(t)') since this will
    # default to assuming that t is complex.  So if the symbol has the
    # same representation, convert to the desired one.
    var = sym.Symbol(str(t))
    expr = expr.replace(var, t)

    orig_expr = expr

    # sym.rewrite(sym.exp) can create exp(log...)
    if expr.has(sym.sin):
        expr = expr.replace(lambda expr: expr.is_Function and expr.func == sym.sin,
                            lambda expr: expr.rewrite(sym.exp))
    if expr.has(sym.cos):
        expr = expr.replace(lambda expr: expr.is_Function and expr.func == sym.cos,
                            lambda expr: expr.rewrite(sym.exp))        

    terms = expr.expand().as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += fourier_term(symsimplify(term), t, f, inverse=inverse)
    except ValueError:
        raise ValueError('Could not compute Fourier transform for ' + str(orig_expr))

    fourier_cache[key] = result
    return result


def inverse_fourier_transform(expr, f, t, evaluate=True):
    """Compute bilateral inverse Fourier transform of expr.

    Undefined functions such as V(f) are converted to v(t)

    This also handles some expressions that do not really have an
    inverse Fourier transform, such as a, cos(a * f), sin(a * f), exp(I *
    a * f)."""
    
    if expr.has(t):
        raise ValueError('Cannot inverse Fourier transform for expression %s that depends on %s' % (expr, t))
    
    result = fourier_transform(expr, t, f, inverse=True, evaluate=evaluate)
    return symsimplify(result)


def FT(expr, t, f, **assumptions):
    """Compute bilateral Fourier transform of expr.

    Undefined functions such as v(t) are converted to V(f)

    This also handles some expressions that do not really have a Fourier
    transform, such as a, cos(a * t), sin(a * t), exp(I * a * t)."""    
    
    return fourier_transform(expr, t, f, **assumptions)


def IFT(expr, f, t, **assumptions):
    """Compute bilateral inverse Fourier transform of expr.

    Undefined functions such as V(f) are converted to v(t)

    This also handles some expressions that do not really have an
    inverse Fourier transform, such as a, cos(a * f), sin(a * f), exp(I *
    a * f)."""    
    
    return inverse_fourier_transform(expr, f, t, **assumptions)

"""This module provides support for Laplace transforms.  It acts as a
quantity for SymPy's Laplace transform.  It calculates the unilateral
Laplace transform using:

F(s) = lim_{t_0\rightarrow 0} \int_{-t_0}^{\infty} f(t) e^{-s t} dt

In comparison, SymPy uses:

F(s) = \int_{0}^{\infty} f(t) e^{-s t} dt

The latter gives 0.5 for the Laplace transform of DiracDelta(t)
whereas the former version gives 1.  Note, SymPy is inconsistent in
that it gives DiracDelta(t) for the inverse Laplace transform of 1.

Another difference with this implementation is that it will transform
undefined functions such as v(t) to V(s).

These functions are for internal use by Lcapy.  

Copyright 2016--2021 Michael Hayes, UCECE

"""

from .ratfun import Ratfun
from .sym import sympify, simplify, AppliedUndef
from .utils import factor_const, scale_shift, as_sum_terms
import sympy as sym

__all__ = ('LT', 'laplace_transform')

laplace_cache = {}


def laplace_limits(expr, t, s, tmin, tmax):
    
    F = sym.integrate(expr * sym.exp(-s * t), (t, tmin, tmax))

    if not F.has(sym.Integral):
        return F

    if not F.is_Piecewise:
        raise ValueError('Could not compute Laplace transform for ' + str(expr))

    F, cond = F.args[0]
    if F.has(sym.Integral):
        raise ValueError('Could not compute Laplace transform for ' + str(expr))

    return F


def laplace_0minus(expr, t, s):
    
    t0 = sym.symbols('t0', negative=True, real=True)

    F = laplace_limits(expr, t, s, t0, sym.oo)
    return sym.limit(F, t0, 0)


def laplace_0(expr, t, s):

    return laplace_limits(expr, t, s, 0, sym.oo)


def laplace_func(expr, t, s, inverse=False):

    if not isinstance(expr, AppliedUndef):
        raise ValueError('Expecting function for %s' % expr)

    scale, shift = scale_shift(expr.args[0], t)    

    ssym = sympify(str(s))
    
    # Convert v(t) to V(s), etc.
    name = expr.func.__name__
    if inverse:
        func = name[0].lower() + name[1:] + '(%s)' % s
    else:
        func = name[0].upper() + name[1:] + '(%s)' % s

    result = sympify(func).subs(ssym, s / scale) / abs(scale)

    if shift != 0:
        result = result * sym.exp(s * shift / scale)    
    return result


def laplace_integral(expr, t, s):

    const, expr = factor_const(expr, t)

    if len(expr.args) != 2:
        raise ValueError('Cannot compute Laplace transform of %s' % expr)

    integrand = expr.args[0]
    
    if not isinstance(expr, sym.Integral):
        raise ValueError('Cannot compute Laplace transform of %s' % expr)

    if len(expr.args[1]) != 3:
        raise ValueError('Require definite integral')
    
    var = expr.args[1][0]
    limits = expr.args[1][1:]
    const2, expr2 = factor_const(integrand, var)

    if (expr2.is_Function and
        expr2.args[0] == t - var and limits[0] == 0 and limits[1] == sym.oo):
        return const2 * laplace_term(expr2.subs(t - var, t), t, s) / s

    # Look for convolution integral
    # TODO, handle convolution with causal functions.
    if (limits[0] != -sym.oo) or (limits[1] != sym.oo):
        raise ValueError('Need indefinite limits for %s' % expr)
   
    if ((len(expr.args) != 2) or not expr2.is_Mul or
        not expr2.args[0].is_Function or not expr2.args[1].is_Function):
        raise ValueError('Need integral of product of two functions: %s' % expr)

    f1 = expr2.args[0]
    f2 = expr2.args[1]    
    # TODO: apply similarity theorem if have f(a * tau) etc.

    if (f1.args[0] == var and f2.args[0] == t - var):
        F1 = laplace_term(f1, var, s)
        F2 = laplace_term(f2.subs(t - var, t), t, s)
    elif (f2.args[0] == var and f1.args[0] == t - var):
        F1 = laplace_term(f1.subs(t - var, t), t, s)
        F2 = laplace_term(f2, var, s)
    else:            
        raise ValueError('Cannot recognise convolution: %s' % expr)
    
    return const2 * F1 * F2


def laplace_derivative_undef(expr, t, s):
    
    if not isinstance(expr, sym.Derivative):
        raise ValueError('Cannot compute Laplace transform of %s' % expr)
    
    if (not isinstance(expr.args[0], AppliedUndef) and
        expr.args[1][0] != t):
        raise ValueError('Cannot compute Laplace transform of %s' % expr)

    ssym = sympify(str(s))    
    name = expr.args[0].func.__name__    
    func1 = name[0].upper() + name[1:] + '(%s)' % str(ssym)    
    return sympify(func1).subs(ssym, s) * s ** expr.args[1][1]


def laplace_sin_cos(expr, t, s):

    # Sympy sometimes has problems with this...
    
    if not expr.is_Function or expr.func not in (sym.sin, sym.cos):
        raise ValueError('Expression not sin or cos')
    arg = expr.args[0]
    a, b = scale_shift(arg, t)

    d = s**2 + a**2

    if expr.func is sym.sin:
        return (a * sym.cos(b) + s * sym.sin(b)) / d
    else:
        return (-a * sym.sin(b) + s * sym.cos(b)) / d        

def laplace_term(expr, t, s):

    # Unilateral LT ignores expr for t < 0 so remove Piecewise.
    if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
        expr = expr.args[0].args[0]
    
    const, expr = factor_const(expr, t)

    terms = expr.as_ordered_terms()
    if len(terms) > 1:
        result = 0
        for term in terms:
            result += laplace_term(term, t, s)
        return const * result

    tsym = sympify(str(t))
    expr = expr.replace(tsym, t)

    if expr.has(sym.Integral):
        return laplace_integral(expr, t, s) * const

    if expr.is_Function and expr.func in (sym.sin, sym.cos) and expr.args[0].has(t):
        return laplace_sin_cos(expr, t, s) * const
    
    if expr.has(AppliedUndef):

        if expr.has(sym.Derivative):
            return laplace_derivative_undef(expr, t, s) * const    

        factors = expr.as_ordered_factors()
        if len(factors) == 1:
            return const * laplace_func(factors[0], t, s)
        elif len(factors) > 2:
            raise ValueError('TODO: cannot handle product %s' % expr)

        foo = factors[1]
        if foo.is_Function and foo.func == sym.exp and foo.args[0].has(t):
            scale, shift = scale_shift(foo.args[0], t)
            if shift == 0: 
                result = laplace_func(factors[0], t, s)
                return const * result.subs(s, s - scale)
        raise ValueError('TODO: cannot handle product %s' % expr)

    if expr.has(sym.Heaviside(t)):
        return laplace_0(expr.replace(sym.Heaviside(t), 1), t, s) * const

    if expr.has(sym.DiracDelta) or expr.has(sym.Heaviside):
        try:
            return laplace_0minus(expr, t, s) * const
        except ValueError:
            pass

    return laplace_0(expr, t, s) * const


def laplace_transform(expr, t, s, evaluate=True):
    """Compute unilateral Laplace transform of expr with lower limit 0-.

    Undefined functions such as v(t) are converted to V(s)."""

    if expr.is_Equality:
        return sym.Eq(laplace_transform(expr.args[0], t, s),
                      laplace_transform(expr.args[1], t, s))

    if not evaluate:
        t0 = sympify('t0')
        result = sym.Limit(sym.Integral(expr * sym.exp(-s * t), (t, t0, sym.oo)),
                           t0, 0, dir='-')
        return result

    # Unilateral LT ignores expr for t < 0 so remove Piecewise.
    if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
        expr = expr.args[0].args[0]
    
    const, expr = factor_const(expr, t)    
    
    key = (expr, t, s)
    if key in laplace_cache:
        return const * laplace_cache[key]

    if expr.has(s):
        raise ValueError('Cannot Laplace transform for expression %s that depends on %s' % (expr, s))
    
    # The variable may have been created with different attributes,
    # say when using sympify('Heaviside(t)') since this will
    # default to assuming that t is complex.  So if the symbol has the
    # same representation, convert to the desired one.

    var = sym.Symbol(str(t))
    if isinstance(expr, Expr):
        expr = expr.expr
    else:
        expr = sympify(expr)

    expr = expr.replace(var, t)        
        
    # Unilateral LT ignores expr for t < 0 so remove Piecewise.
    if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
        expr = expr.args[0].args[0]

    # expand can make a mess of expressions with exponentials
    #expr = sym.expand(expr)        
    terms = expr.as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += laplace_term(term, t, s)
    except ValueError:
        raise

    result = result.simplify()
    laplace_cache[key] = result
    return const * result



def inverse_laplace_make(t, const, cresult, uresult, **assumptions):

    result = const * (cresult + uresult)
    
    if assumptions.get('dc', False):
        free_symbols = set([symbol.name for symbol in result.free_symbols])
        if 't' in free_symbols:
            raise ValueError('Something wonky going on, expecting dc.'
                             ' Perhaps have capacitors in series?')
    
    elif assumptions.get('ac', False):

        if cresult != 0:
            # Quietly drop DiracDelta term.  This will appear when
            # trying to determine the current through a capacitor
            # when a cosinusoidal voltage source is applied.
            if cresult.has(sym.DiracDelta):
                result = uresult
            else:
                raise ValueError('Inverse laplace transform weirdness for %s'
                                 ' with is_ac True' % result)
            
        # TODO, perform more checking of the result.
        
    elif not assumptions.get('causal', False):

        if uresult != 0:
            # Cannot determine result for t < 0
            result = sym.Piecewise((const * (uresult + cresult), t >= 0))
        else:
            result = const * cresult
        
    return result


def LT(expr, t, s, **assumptions):
    """Compute unilateral Laplace transform of expr with lower limit 0-.

    Undefined functions such as v(t) are converted to V(s)."""
    
    return laplace_transform(expr, t, s, **assumptions)


from .expr import Expr

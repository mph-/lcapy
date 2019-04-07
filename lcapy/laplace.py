"""This module provides support for Laplace transforms.  It acts as a
wrapper for SymPy's Laplace transform.  It calculates the unilateral
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

Copyright 2016--2019 Michael Hayes, UCECE

"""

from .ratfun import Ratfun
from .sym import sympify
from .utils import factor_const, scale_shift
import sympy as sym

laplace_cache = {}
inverse_laplace_cache = {}


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

    if not isinstance(expr, sym.function.AppliedUndef):
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

    if len(expr.args) != 1:
        raise ValueError('Cannot compute Laplace transform of %s' % expr)

    expr = expr.args[0]
    
    if not isinstance(expr, sym.Integral):
        raise ValueError('Cannot compute Laplace transform of %s' % expr)

    # Look for convolution integral
    
    var = expr.args[1][0]
    if (expr.args[1][1] != -sym.oo) or (expr.args[1][2] != sym.oo):
        raise ValueError('Need indefinite limits for %s' % expr)
    
    const2, expr = factor_const(expr.args[0], t)
    if ((len(expr.args) != 2)
        or (not isinstance(expr.args[0], sym.function.AppliedUndef))
        or (not isinstance(expr.args[1], sym.function.AppliedUndef))):
        raise ValueError('Need integral of two functions: %s' % expr)        

    f1 = expr.args[0]
    f2 = expr.args[1]    

    # TODO: apply similarity theorem if have f(a tau) etc.
    
    if ((f1.args[0] != var or f2.args[0] != t - var)
        and (f2.args[0] != var or f1.args[0] != t - var)):
        raise ValueError('Cannot recognise convolution: %s' % expr)

    ssym = sympify(str(s))
    
    name = f1.func.__name__
    func1 = name[0].upper() + name[1:] + '(%s)' % str(ssym)

    name = f2.func.__name__
    func2 = name[0].upper() + name[1:] + '(%s)' % str(ssym)    

    F1 = sympify(func1).subs(ssym, s)
    F2 = sympify(func2).subs(ssym, s)
    
    return F1 * F2

def laplace_term(expr, t, s):

    const, expr = factor_const(expr, t)

    tsym = sympify(str(t))
    expr = expr.replace(tsym, t)

    if expr.has(sym.Integral):
        return laplace_integral(expr, t, s) * const
    
    if expr.has(sym.function.AppliedUndef):

        rest = sym.S.One
        for factor in expr.as_ordered_factors():
            if isinstance(factor, sym.function.AppliedUndef):
                result = laplace_func(factor, t, s)
            else:
                if factor.has(t):
                    raise ValueError('TODO: need derivative of undefined'
                                     ' function for %s' % factor)
                rest *= factor
        return result * rest * const

    if expr.has(sym.Heaviside(t)):
        return laplace_0(expr.replace(sym.Heaviside(t), 1), t, s) * const

    if expr.has(sym.DiracDelta) or expr.has(sym.Heaviside):
        try:
            return laplace_0minus(expr, t, s) * const
        except ValueError:
            pass

    return laplace_0(expr, t, s) * const


def laplace_transform(expr, t, s):
    """Compute unilateral Laplace transform of expr with lower limit 0-.

    Undefined functions such as v(t) are converted to V(s)

    """

    key = (expr, t, s)
    if key in laplace_cache:
        return laplace_cache[key]

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

    # Unilateral LT ignores expr for t < 0 so
    # but barfs on a Piecewise so handle case here.
    expr = expr.replace(var, t)        
    if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
        expr = expr.args[0].args[0]

    terms = expr.as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += laplace_term(term, t, s)
    except ValueError:
        raise

    result = result.simplify()
    laplace_cache[key] = result
    return result


def inverse_laplace_ratfun(expr, s, t):

    N, D, delay = Ratfun(expr, s).as_ratfun_delay()
    # The delay should be zero

    Q, M = sym.div(N, D, s)

    result1 = sym.S.Zero

    if Q:
        Qpoly = sym.Poly(Q, s)        
        C = Qpoly.all_coeffs()
        for n, c in enumerate(C):
            result1 += c * sym.diff(sym.DiracDelta(t), t, len(C) - n - 1)

    expr = M / D
    for factor in expr.as_ordered_factors():
        if factor == sym.oo:
            return factor

    sexpr = Ratfun(expr, s)
    P = sexpr.poles()
    result2 = sym.S.Zero

    P2 = P.copy()

    for p in P2:

        # Number of occurrences of the pole.
        N = P2[p]

        if N == 0:
            continue

        f = s - p

        if N == 1:
            r = sexpr.residue(p, P)

            pc = p.conjugate()
            if pc != p and pc in P:
                # Remove conjugate from poles and process pole with its
                # conjugate.  Unfortunately, for symbolic expressions
                # we cannot tell if a quadratic has two real poles,
                # a repeat real pole, or a complex conjugate pair of poles.
                P2[pc] = 0
                
                p_re = sym.re(p)
                p_im = sym.im(p)
                r_re = sym.re(r)
                r_im = sym.im(r)
                et = sym.exp(p_re * t)
                result2 += 2 * r_re * et * sym.cos(p_im * t)
                result2 -= 2 * r_im * et * sym.sin(p_im * t)
            else:
                result2 += r * sym.exp(p * t)
            continue

        # Handle repeated poles.
        expr2 = expr * f ** N
        for n in range(1, N + 1):
            m = N - n
            r = sym.limit(
                sym.diff(expr2, s, m), s, p) / sym.factorial(m)
            result2 += r * sym.exp(p * t) * t**(n - 1)

    # result1 is a sum of Dirac deltas and its derivatives so is known
    # to be causal.

    return result1, result2


def dummyvar(intnum=0):
    if intnum == 0:
        return sympify('tau')
    else:
        return sympify('tau_%d' % intnum)    


def inverse_laplace_product(expr, s, t, **assumptions):

    # Handle expressions with a function of s, e.g., V(s) * Y(s), V(s)
    # / s etc.

    if assumptions.get('causal', False):
        # Assume that all functions are causal in the expression.
        t1 = sym.S.Zero
        t2 = t
    else:
        t1 = -sym.oo
        t2 = sym.oo        
    
    const, expr = factor_const(expr, s)

    factors = expr.as_ordered_factors()
    if len(factors) < 2:
        raise ValueError('Expression does not have multiple factors: %s' % expr)

    if isinstance(factors[1], sym.function.AppliedUndef):
        # Try to expose more simple cases, e.g. (R + s * L) * V(s)
        terms = factors[0].as_ordered_terms()
        if len(terms) >= 2:
            result = sym.S.Zero
            for term in terms:
                result += inverse_laplace_product(factors[1] * term, s, t)
            return result * const

    result1, result2 = inverse_laplace_term1(factors[0], s, t)
    result = result1 + result2

    intnum = 0
    for m in range(len(factors) - 1):
        if m == 0 and isinstance(factors[1], sym.function.AppliedUndef):
            # Note, as_ordered_factors puts powers of s before the functions.
            if factors[0] == s:
                # Handle differentiation
                # Convert s * V(s) to d v(t) / dt                        
                result = laplace_func(factors[1], s, t, True)            
                result = sym.Derivative(result, t)
                continue
            elif factors[0].is_Pow and factors[0].args[0] == s and factors[0].args[1] > 0:
                # Handle higher order differentiation
                # Convert s ** 2 * V(s) to d^2 v(t) / dt^2
                result = laplace_func(factors[1], s, t, True)            
                result = sym.Derivative(result, t, factors[0].args[1])
                continue                
            elif factors[0].is_Pow and factors[0].args[0] == s and factors[0].args[1] == -1:
                # Handle integration  1 / s * V(s)
                tau = dummyvar(intnum)
                intnum += 1
                result = laplace_func(factors[1], s, tau, True)
                result = sym.Integral(result, (tau, t1, t))
                continue                
        # Convert product to convolution
        tau = dummyvar(intnum)
        intnum += 1
        result1, result2 = inverse_laplace_term1(factors[m + 1], s, t)
        expr2 = result1 + result2
        result = sym.Integral(result.subs(t, t - tau) * expr2.subs(t, tau),
                              (tau, t1, t2))
    
    return result * const


def delay_factor(expr, var):

    delay = sym.S.Zero    
    rest = sym.S.One
    
    for f in expr.as_ordered_factors():
        b, e = f.as_base_exp()
        if b == sym.E and e.is_polynomial(var):
            p = sym.Poly(e, var)
            c = p.all_coeffs()
            if p.degree() == 1:
                delay -= c[0]
                if c[1] != 0:
                    rest *= sym.exp(c[1])
                continue

        rest *= f
    return rest, delay


def inverse_laplace_sympy(expr, s, t):

    # This barfs when needing to generate Dirac deltas
    from sympy.integrals.transforms import inverse_laplace_transform
    result = inverse_laplace_transform(expr, t, s)
    
    if result.has(sym.InverseLaplaceTransform):
        raise ValueError('Cannot determine inverse Laplace'
                         ' transform of %s with sympy' % expr)
    return result


def inverse_laplace_term1(expr, s, t, **assumptions):

    const, expr = factor_const(expr, s)

    if isinstance(expr, sym.function.AppliedUndef):
        # Handle V(s), 3 * V(s) etc.  If causal is True it is assumed
        # that the unknown functions are causal.  Note laplace_func
        # just changes the name so it works as inverse_laplace_func.
        result = laplace_func(expr, s, t, True)
        return result * const, sym.S.Zero
    
    if expr.has(sym.function.AppliedUndef):
        return const * inverse_laplace_product(expr, s, t,
                                               **assumptions), sym.S.Zero

    try:
        # This is the common case.
        result1, result2 = inverse_laplace_ratfun(expr, s, t)
        return const * result1, const * result2
    except:
        pass

    try:
        return sym.S.Zero, const * inverse_laplace_sympy(expr, s, t, **assumptions)
    except:
        pass

    # As last resort see if can convert to convolutions...
    return sym.S.Zero, const * inverse_laplace_product(expr, s, t)
    
    
def inverse_laplace_term(expr, s, t, **assumptions):

    expr, delay = delay_factor(sym.simplify(expr), s)

    result1, result2 = inverse_laplace_term1(expr, s, t, **assumptions)

    if delay != 0:
        result1 = result1.subs(t, t - delay)
        result2 = result2.subs(t, t - delay)

    # TODO, should check for delay < 0.  If so the causal
    # part is no longer causal.

    if assumptions.get('causal', False):
        result2 = result2 * sym.Heaviside(t - delay)
    
    return result1, result2


def inverse_laplace_by_terms(expr, s, t, **assumptions):

    expr = sym.expand(expr)
    terms = expr.as_ordered_terms()

    result1 = sym.S.Zero
    result2 = sym.S.Zero    

    for term in terms:
        part1, part2 = inverse_laplace_term(term, s, t, **assumptions)
        result1 += part1
        result2 += part2        
    return result1, result2


def inverse_laplace_transform(expr, s, t, **assumptions):
    """Calculate inverse Laplace transform of X(s) and return x(t).

    The unilateral Laplace transform cannot determine x(t) for t < 0
    unless given additional information in the way of assumptions.

    The assumptions are:
    dc -- x(t) = constant so X(s) must have the form constant / s
    causal -- x(t) = 0 for t < 0.
    ac -- x(t) = A cos(a * t) + B * sin(b * t)
    """

    # TODO, simplify
    key = (expr, s, t, assumptions.get('dc', False),
           assumptions.get('ac', False),
           assumptions.get('causal', False))
    
    if key in inverse_laplace_cache:
        return inverse_laplace_cache[key]

    if expr.has(t):
        raise ValueError('Cannot inverse Laplace transform for expression %s that depends on %s' % (expr, t))
    
    if assumptions.get('dc', False):
        result = expr * s
            
        free_symbols = set([symbol.name for symbol in result.free_symbols])
        if 's' in free_symbols:
            raise ValueError('Something wonky going on, expecting dc.'
                             ' Perhaps have capacitors in series?')
        return result

    try:
        result1, result2 = inverse_laplace_term(expr, s, t, **assumptions)
    except:
        result1, result2 = inverse_laplace_by_terms(expr, s, t, **assumptions)
        
    # result1 is known to be causal, result2 is unsure
    result = result1 + result2    

    if assumptions.get('ac', False):
        if result1 != 0:
            raise ValueError('Inverse laplace transform weirdness for %s'
                             ' with is_ac True' % expr)
        # TODO, perform more checking of the result.
        
    elif not assumptions.get('causal', False):
        result = sym.Piecewise((result, t >= 0))
        
    inverse_laplace_cache[key] = result
    return result

from .expr import Expr

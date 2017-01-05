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

Copyright 2016 Michael Hayes, UCECE

"""

from lcapy.ratfun import Ratfun
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

def laplace_term(expr, t, s):

    var = sym.Symbol(str(t))
    expr = expr.replace(var, t)

    if expr.has(sym.function.AppliedUndef) and expr.args[0] == t:
        # TODO, handle things like 3 * v(t), a * v(t), 3 * t * v(t), v(t-T),
        # v(4 * a * t), etc.
        if not isinstance(expr, sym.function.AppliedUndef):
            raise ValueError('Could not compute Laplace transform for ' + str(expr))

        # Convert v(t) to V(s), etc.
        name = expr.func.__name__
        name = name[0].upper() + name[1:] + '(s)'
        return sym.sympify(name)

    if expr.has(sym.Heaviside(t)):
        return laplace_0(expr.replace(sym.Heaviside(t), 1), t, s)

    if expr.has(sym.DiracDelta) or expr.has(sym.Heaviside):
        try:
            return laplace_0minus(expr, t, s)
        except ValueError:
            pass

    return laplace_0(expr, t, s)


def laplace_transform(expr, t, s):
    """Compute unilateral Laplace transform of expr with lower limit 0-.

    Undefined functions such as v(t) are converted to V(s)

    """

    key = (expr, t, s)
    if key in laplace_cache:
        return laplace_cache[key]
    
    # The variable may have been created with different attributes,
    # say when using sym.sympify('Heaviside(t)') since this will
    # default to assuming that t is complex.  So if the symbol has the
    # same representation, convert to the desired one.

    var = sym.Symbol(str(t))
    if hasattr(expr, 'expr'):
        expr = expr.expr
    else:
        expr = sym.sympify(expr)

    # Unilateral LT ignores expr for t < 0 so
    # but barfs on a Piecewise so handle case here.
    expr = expr.replace(var, t)        
    if expr.is_Piecewise and expr.args[0].args[1] == t >= 0:
        expr = expr.args[0].args[0]

    terms = expr.as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += laplace_term(term, t, s)
    except ValueError:
        raise ValueError('Could not compute Laplace transform for ' + str(expr))

    result = result.simplify()
    laplace_cache[key] = result
    return result


def inverse_laplace_ratfun(expr, s, t, **assumptions):

    N, D, delay = Ratfun(expr, s).as_ratfun_delay()

    Q, M = N.div(D)

    result1 = 0

    # Delayed time.
    td = t - delay

    if Q:
        C = Q.all_coeffs()
        for n, c in enumerate(C):
            result1 += c * sym.diff(sym.DiracDelta(td), t, len(C) - n - 1)

    expr = M / D
    for factor in expr.as_ordered_factors():
        if factor == sym.oo:
            return factor

    sexpr = Ratfun(expr, s)
    P = sexpr.poles()
    result2 = 0

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
                etd = sym.exp(p_re * td)
                result2 += 2 * r_re * etd * sym.cos(p_im * td)
                result2 -= 2 * r_im * etd * sym.sin(p_im * td)
            else:
                result2 += r * sym.exp(p * td)
            continue

        # Handle repeated poles.
        expr2 = expr * f ** N
        for n in range(1, N + 1):
            m = N - n
            r = sym.limit(
                sym.diff(expr2, s, m), s, p) / sym.factorial(m)
            result2 += r * sym.exp(p * td) * td**(n - 1)

    if assumptions.get('ac', False):
        return result1 + result2

    if assumptions.get('causal', False):
        return result1 + result2 * sym.Heaviside(td)

    if delay != 0:
        result2 *= sym.Heaviside(td)
                
    return sym.Piecewise((result1 + result2, t >= 0))


def inverse_laplace_special(expr, s, t, **assumptions):

    if expr.has(sym.function.AppliedUndef) and expr.args[0] == s:

        # TODO, handle things like 3 * V(s), a * V(s), 3 * s * V(s), etc.
        if isinstance(expr, sym.function.AppliedUndef):
            # Convert V(s) to v(t), etc.
            name = expr.func.__name__
            name = name[0].lower() + name[1:] + '(t)'
            return sym.sympify(name)
    
    raise ValueError('Could not compute inverse Laplace transform for ' + str(expr))
    
    
def inverse_laplace_term(expr, s, t, **assumptions):

    try:
        return inverse_laplace_ratfun(expr, s, t, **assumptions)
    except:
        pass

    try:
        return inverse_laplace_special(expr, s, t, **assumptions)
    except:
        pass
    
    try:
        # Try splitting into partial fractions to help sympy.
        expr = Ratfun(expr, s).partfrac()
    except:
        pass
            
    # This barfs when needing to generate Dirac deltas
    from sympy.integrals.transforms import inverse_laplace_transform
    result = inverse_laplace_transform(expr, t, s)
    
    if result.has(sym.InverseLaplaceTransform):
        raise ValueError('Cannot determine inverse Laplace'
                         ' transform of %s with sympy' % expr)
   
    return result

def inverse_laplace_transform(expr, s, t, **assumptions):

    # TODO, simplify
    key = (expr, s, t, assumptions.get('dc', False),
           assumptions.get('ac', False),
           assumptions.get('causal', False))
    
    if key in inverse_laplace_cache:
        return inverse_laplace_cache[key]

    if assumptions.get('dc', False):
        result = expr * s
            
        free_symbols = set([symbol.name for symbol in result.free_symbols])
        if 's' in free_symbols:
            raise ValueError('Something wonky going on, expecting dc.'
                             ' Perhaps have capacitors in series?')
        return result

    try:
        result = inverse_laplace_term(expr, s, t, **assumptions)
    except:
        terms = expr.as_ordered_terms()
        result = 0

        for term in terms:
            result += inverse_laplace_term(term, s, t)

        # Perhaps could cache terms

    if not assumptions.get('causal', False):
        result = sym.Piecewise((result, t >= 0))
        
    inverse_laplace_cache[key] = result
    return result

    

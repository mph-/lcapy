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


def inverse_laplace_ratfun(expr, s, t):

    N, D, delay = Ratfun(expr, s).as_ratfun_delay()
    # The delay should be zero

    Q, M = N.div(D)

    result1 = sym.sympify(0)

    if Q:
        C = Q.all_coeffs()
        for n, c in enumerate(C):
            result1 += c * sym.diff(sym.DiracDelta(t), t, len(C) - n - 1)

    expr = M / D
    for factor in expr.as_ordered_factors():
        if factor == sym.oo:
            return factor

    sexpr = Ratfun(expr, s)
    P = sexpr.poles()
    result2 = sym.sympify(0)

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


def inverse_laplace_special(expr, s, t):

    if not expr.has(sym.function.AppliedUndef):
        raise ValueError('Could not compute inverse Laplace transform for ' + str(expr))

    ssym = sym.sympify(str(s))
    expr = expr.subs(ssym, s)

    rest = sym.sympify(1)
    for factor in expr.as_ordered_factors():
        if isinstance(factor, sym.function.AppliedUndef):
            if factor.args[0] != s:
                raise ValueError('Weird function %s not of s' % factor)
                
            # Convert V(s) to v(t), etc.
            name = factor.func.__name__
            tsym = sym.sympify(str(t))
            func = name[0].lower() + name[1:] + '(%s)' % str(tsym)
            result = sym.sympify(func).subs(tsym, t)
        else:
            if factor.has(s):
                raise ValueError('TODO: need derivative of undefined'
                                 ' function for %s' % factor)
            rest *= factor

    return rest * result


def delay_factor(expr, var):

    delay = sym.sympify(0)    
    rest = sym.sympify(1)
    
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

    try:
        # Handle V(s), etc.
        return sym.sympify(0), inverse_laplace_special(expr, s, t)
    except:
        pass

    try:
        # This is the common case.
        return inverse_laplace_ratfun(expr, s, t)
    except:
        pass

    try:
        return sym.sympify(0), inverse_laplace_sympy(expr, s, t)
    except:
        pass    
    
    raise ValueError('Cannot determine inverse Laplace'
                     ' transform of %s with sympy' % expr)    


def inverse_laplace_term(expr, s, t, **assumptions):

    expr, delay = delay_factor(expr, s)

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

    terms = expr.as_ordered_terms()

    result1 = sym.sympify(0)
    result2 = sym.sympify(0)    

    for term in terms:
        part1, part2 = inverse_laplace_term(term, s, t, **assumptions)
        result1 += part1
        result2 += part2        
    return result1, result2


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
        result1, result2 = inverse_laplace_term(expr, s, t, **assumptions)
    except:
        result1, result2 = inverse_laplace_by_terms(expr, s, t, **assumptions)
        
    # result1 is known to be causal, result2 is unsure

    if assumptions.get('ac', False):
        if result1 != 0:
            raise ValueError('Inverse laplace transform weirdness for %s'
                             ' with is_ac True' % expr)
        result = result1 + result2
    elif not assumptions.get('causal', False):
        result = sym.Piecewise((result1 + result2, t >= 0))
    else:
        result = result1 + result2        
        
    inverse_laplace_cache[key] = result
    return result

    

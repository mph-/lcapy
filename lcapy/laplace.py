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

__all__ = ('LT', 'ILT')

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


def inverse_laplace_damped_sin(expr, s, t, **assumptions):

    ncoeffs, dcoeffs = expr.coeffs()
    K = ncoeffs[0] / dcoeffs[0]

    ncoeffs = [(c / ncoeffs[0]) for c in ncoeffs]
    dcoeffs = [(c / dcoeffs[0]) for c in dcoeffs]        

    if len(ncoeffs) > 3 or len(dcoeffs) > 3:
        raise ValueError('Not a second-order response')
    
    omega0 = sym.sqrt(dcoeffs[2])
    zeta = dcoeffs[1] / (2 * omega0)

    if zeta.is_constant() and zeta > 1:
        print('Warning: expression is overdamped')

    sigma1 = (zeta * omega0).simplify()
    omega1 = (omega0 * sym.sqrt(1 - zeta**2)).simplify()
    K = (K / omega1).simplify()

    E = sym.exp(-sigma1 * t)
    S = sym.sin(omega1 * t)

    h = K * E * S

    # If overdamped
    #h = K * sym.exp(-sigma1 * t) * sym.sinh(omega0 * mu * t)
        
    if len(ncoeffs) == 1:
        return sym.S.Zero, h

    C = sym.cos(omega1 * t)
    kCd = omega1
    kSd = -sigma1
    hd = K * E * (kCd * C + kSd * S)
    
    if len(ncoeffs) == 2:
        return sym.S.Zero, K * E * (kCd * C + (ncoeffs[1] + kSd) * S)

    kCdd = -2 * omega1 * sigma1
    kSdd = sigma1**2 - omega1**2

    G = K * E * ((kCdd + ncoeffs[1] * kCd) * C + (kSdd + ncoeffs[1] * kSd + ncoeffs[2]) * S)
    
    return K * kCd * sym.DiracDelta(t), G


def inverse_laplace_ratfun(expr, s, t, **assumptions):

    sexpr = Ratfun(expr, s)

    damping = assumptions.get('damping', None)

    if assumptions.get('damped_sin', False):
        if sexpr.degree == 2:
            return inverse_laplace_damped_sin(sexpr, s, t, **assumptions)
        #if False and sexpr.degree == 3 and Ratfun(expr * s).degree == 2:
        #    return inverse_laplace_damped_sin3(sexpr, s, t, **assumptions)

    Q, M, D, delay, undef = sexpr.as_QMD()

    cresult = sym.S.Zero

    if Q:
        Qpoly = sym.Poly(Q, s)        
        C = Qpoly.all_coeffs()
        for n, c in enumerate(C):
            cresult += c * sym.diff(sym.DiracDelta(t), t, len(C) - n - 1)

    expr = M / D
    for factor in expr.as_ordered_factors():
        if factor == sym.oo:
            return factor

    sexpr = Ratfun(expr, s)
    poles = sexpr.poles(damping=damping)
    polesdict = {}
    for pole in poles:
        polesdict[pole.expr] = pole.n
    
    uresult = sym.S.Zero

    for pole in poles:

        p = pole.expr

        # Number of occurrences of the pole.
        o = polesdict[p]        

        if o == 0:
            continue

        if o == 1:
            pc = pole.conjugate
            r = sexpr.residue(p, poles)
            
            if pc != p and pc in polesdict:
                # Remove conjugate from poles and process pole with its
                # conjugate.  Unfortunately, for symbolic expressions
                # we cannot tell if a quadratic has two real poles,
                # a repeated real pole, or a complex conjugate pair of poles.

                polesdict[pc] -= 1
                
                p_re = sym.re(p)
                p_im = sym.im(p)
                r_re = sym.re(r)
                r_im = sym.im(r)
                et = sym.exp(p_re * t)
                uresult += 2 * r_re * et * sym.cos(p_im * t)
                uresult -= 2 * r_im * et * sym.sin(p_im * t)
            else:
                uresult += r * sym.exp(p * t)
            continue

        # Handle repeated poles.
        expr2 = expr * (s - p) ** o
        expr2 = expr2.simplify()        
        for n in range(1, o + 1):
            m = o - n
            r = sym.limit(
                sym.diff(expr2, s, m), s, p) / sym.factorial(m)
            uresult += r * sym.exp(p * t) * t**(n - 1)

    # cresult is a sum of Dirac deltas and its derivatives so is known
    # to be causal.

    return cresult, uresult


def dummyvar(intnum=0, var='tau'):
    if intnum == 0:
        return sympify(var, real=True)
    else:
        return sympify('%s_%d' % (var, intnum), real=True)    


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

    if (len(factors) > 2 and not
        # Help s * 1 / (s + R * C) * I(s)
        isinstance(factors[1], AppliedUndef) and
        isinstance(factors[2], AppliedUndef)):
        factors = [factors[0], factors[2], factors[1]] + factors[3:]
    
    if isinstance(factors[1], AppliedUndef):
        # Try to expose more simple cases, e.g. (R + s * L) * V(s)
        terms = factors[0].as_ordered_terms()
        if len(terms) >= 2:
            result = sym.S.Zero
            for term in terms:
                result += inverse_laplace_product(factors[1] * term, s, t)
            return result * const

    cresult, uresult = inverse_laplace_term1(factors[0], s, t)
    result = cresult + uresult

    intnum = 0
    for m in range(len(factors) - 1):
        if m == 0 and isinstance(factors[1], AppliedUndef):
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
                tau = dummyvar(intnum, 'tau')
                intnum += 1
                result = laplace_func(factors[1], s, tau, True)
                result = sym.Integral(result, (tau, t1, t))
                continue                
        # Convert product to convolution
        tau = dummyvar(intnum, 'tau')
        intnum += 1
        cresult, uresult = inverse_laplace_term1(factors[m + 1], s, t)
        expr2 = cresult + uresult
        result = sym.Integral(result.subs(t, t - tau) * expr2.subs(t, tau),
                              (tau, t1, t2))
    
    return result * const


def inverse_laplace_power(expr, s, t, **assumptions):

    # Handle expressions with a power of s.
    if not (expr.is_Pow and expr.args[0] == s):
        raise ValueError('Expression %s is not a power of s' % expr)
    exponent = expr.args[1]

    # Have many possible forms; the common ones are:
    # s**a, s**-a, s**(1+a), s**(1-a), s**-(1+a), s**(a-1)
    # Cannot tell if 1-a is positive.

    if exponent.is_positive:
        # Unfortunately, SymPy does not seem to support fractional
        # derivatives...
        return sym.Derivative(sym.DiracDelta(t), t, exponent, evaluate=False)

    if exponent.is_negative:
        return sym.Pow(t, -exponent - 1) / sym.Gamma(-exponent)

    raise ValueError('Cannot determine sign of exponent for %s' % expr)

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

    if isinstance(expr, AppliedUndef):
        # Handle V(s), 3 * V(s) etc.  If causal is True it is assumed
        # that the unknown functions are causal.  Note laplace_func
        # just changes the name so it works as inverse_laplace_func.
        result = laplace_func(expr, s, t, True)
        return result * const, sym.S.Zero
    
    if expr.has(AppliedUndef):
        return const * inverse_laplace_product(expr, s, t,
                                               **assumptions), sym.S.Zero

    try:
        # This is the common case.
        cresult, uresult = inverse_laplace_ratfun(expr, s, t, **assumptions)
        return const * cresult, const * uresult
    except:
        pass

    try:
        return sym.S.Zero, const * inverse_laplace_sympy(expr, s, t, **assumptions)
    except:
        pass

    if expr.is_Pow and expr.args[0] == s:
        return sym.S.Zero, const * inverse_laplace_power(expr, s, t)
    
    # As last resort see if can convert to convolutions...
    return sym.S.Zero, const * inverse_laplace_product(expr, s, t)
    
    
def inverse_laplace_term(expr, s, t, **assumptions):

    expr, delay = delay_factor(expr, s)

    cresult, uresult = inverse_laplace_term1(expr, s, t, **assumptions)

    if delay != 0:
        cresult = cresult.subs(t, t - delay)
        uresult = uresult.subs(t, t - delay)

    # TODO, should check for delay < 0.  If so the causal
    # part is no longer causal.

    if assumptions.get('causal', False):
        uresult = uresult * sym.Heaviside(t - delay)
    
    return cresult, uresult


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


def inverse_laplace_transform1(expr, s, t, verbatim=False, cache_lookup=True,
                               **assumptions):
    """If verbatim is True, do not rewrite expression to assist
    inverse Laplace transform evaluation.

    If cache_lookup is True, try looking for previously cached result.

    """
    
    const, expr = factor_const(expr, s)
    
    key = (expr, s, t,
           assumptions.get('causal', False),                      
           assumptions.get('damping', None),           
           assumptions.get('damped_sin', None))
    
    if key in inverse_laplace_cache:
        cresult, uresult = inverse_laplace_cache[key]        
        return const, cresult, uresult

    if verbatim:
        terms = expr.as_ordered_terms()
    else:
        # This will break how the user constructs the expression but
        # it is more likely to produce a tidy result.
        terms = as_sum_terms(expr, s)
    if len(terms) == 1:
        cresult, uresult = inverse_laplace_term(terms[0], s, t, **assumptions)
    else:
        cresult = sym.S.Zero
        uresult = sym.S.Zero    

        for term in terms:
            const1, part1, part2 = inverse_laplace_transform1(term, s, t, **assumptions)
            cresult += const1 * part1
            uresult += const1 * part2

    # cresult is known to be causal, uresult is unsure            
    inverse_laplace_cache[key] = cresult, uresult
    return const, cresult, uresult


def inverse_laplace_transform(expr, s, t, **assumptions):
    """Calculate inverse Laplace transform of X(s) and return x(t).

    The unilateral Laplace transform cannot determine x(t) for t < 0
    unless given additional information in the way of assumptions.

    The assumptions are:
    dc -- x(t) = constant so X(s) must have the form constant / s
    causal -- x(t) = 0 for t < 0.
    ac -- x(t) = A cos(a * t) + B * sin(b * t)
    """

    if expr.is_Equality:
        return sym.Eq(inverse_laplace_transform(expr.args[0], s, t,
                                                **assumptions),
                      inverse_laplace_transform(expr.args[1], s, t,
                                                **assumptions))    

    if expr.has(t):
        raise ValueError('Cannot inverse Laplace transform for expression %s that depends on %s' % (expr, t))

    const, cresult, uresult = inverse_laplace_transform1(expr, s, t, **assumptions)
    return inverse_laplace_make(t, const, cresult, uresult, **assumptions)    


def LT(expr, t, s, **assumptions):
    """Compute unilateral Laplace transform of expr with lower limit 0-.

    Undefined functions such as v(t) are converted to V(s)."""
    
    return laplace_transform(expr, t, s, **assumptions)


def ILT(expr, s, t, **assumptions):
    """Calculate inverse Laplace transform of X(s) and return x(t).

    The unilateral Laplace transform cannot determine x(t) for t < 0
    unless given additional information in the way of assumptions.

    The assumptions are:
    dc -- x(t) = constant so X(s) must have the form constant / s
    causal -- x(t) = 0 for t < 0.
    ac -- x(t) = A cos(a * t) + B * sin(b * t)
    """
    
    return inverse_laplace_transform(expr, s, t, **assumptions)

from .expr import Expr

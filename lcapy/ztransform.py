"""This module provides support for z transforms.

Copyright 2020 Michael Hayes, UCECE

"""

from .ratfun import Ratfun
from .sym import sympify, simplify, symsymbol
from .utils import factor_const, scale_shift
from .functions import unitimpulse
import sympy as sym

ztransform_cache = {}
inverse_ztransform_cache = {}

# TODO handle convolution -> product


def ztransform_func(expr, n, z, inverse=False):

    if not isinstance(expr, sym.function.AppliedUndef):
        raise ValueError('Expecting function for %s' % expr)

    scale, shift = scale_shift(expr.args[0], n)    

    zsym = sympify(str(z))
    
    # Convert v[n] to V(z), etc.
    name = expr.func.__name__
    if inverse:
        func = name[0].lower() + name[1:] + '(%s)' % z
    else:
        func = name[0].upper() + name[1:] + '(%s)' % z

    if not scale.is_constant():
        raise ValueError('Cannot determine if time-expansion or decimation')

    if scale == 1:
        result = sympify(func).subs(zsym, z)

        if shift != 0:
            result = result * z ** shift
        return result        
    
    if scale.is_integer:
        raise ValueError('Cannot do decimation yet')

    if not scale.is_rational:
        raise ValueError('Cannot handle arbitrary scaling')

    if scale.p != 1:
        raise ValueError('Cannot handle non-integer time-expansion')        
    
    result = sympify(func).subs(zsym, z ** scale.q)

    if shift != 0:
        result = result * z ** shift
    return result


def ztransform_term(expr, n, z):

    const, expr = factor_const(expr, n)

    nsym = sympify(str(n))
    expr = expr.replace(nsym, n)

    if expr.has(sym.function.AppliedUndef):

        rest = sym.S.One
        expr = expr.cancel()
        for factor in expr.as_ordered_factors():
            if isinstance(factor, sym.function.AppliedUndef):
                result = ztransform_func(factor, n, z)
            else:
                if factor.has(n):
                    raise ValueError('TODO: need derivative of undefined'
                                     ' function for %s' % factor)
                rest *= factor
        return result * rest * const

    invz = z ** -1

    result = None
    args = expr.args
    
    if expr.is_Function and expr.func == sym.KroneckerDelta:
        if (args[0] == 0 and args[1] is n) or (args[0] is n and args[1] == 0):
            result = 1
        delay = -(args[1] - args[0] - n)
        if not delay.has(n):
            result = invz ** delay

    elif expr == 1:
        # Unilateral ZT
        result = 1 / (1 - invz)
        
    elif expr.is_Function and expr.func == sym.Heaviside:
        if args[0] is n:
            result = 1 / (1 - invz)
        else:
            delay = args[0] - n
            if not delay.has(n):
                result = invz ** delay * 1 / (1 - invz)

    elif expr.is_Function and expr.func == sym.cos:
        aconst, aexpr = factor_const(args[0], n)
        if aexpr == n:
            result = (1 - sym.cos(aconst) * invz) / (1 - 2 * sym.cos(aconst) * invz + invz ** 2)

    elif expr.is_Function and expr.func == sym.sin:
        aconst, aexpr = factor_const(args[0], n)
        if aexpr == n:
            result = (sym.sin(aconst) * invz) / (1 - 2 * sym.cos(aconst) * invz + invz ** 2)            

    if result is None:
        raise ValueError('Do not know z-transform of %s' % expr)
        
    return result * const


def ztransform(expr, n, z):
    """Compute unilateral z-transform of expr with lower limit 0-.

    Undefined functions such as v[n] are converted to V(z)

    """

    key = (expr, n, z)
    if key in ztransform_cache:
        return ztransform_cache[key]

    if expr.has(z):
        raise ValueError('Cannot Z transform expression %s that depends on %s' % (expr, z))
    
    # The variable may have been created with different attributes,
    # say when using sympify('Heaviside(n)') since this will
    # default to assuming that n is complex.  So if the symbol has the
    # same representation, convert to the desired one.

    var = sym.Symbol(str(n))
    if isinstance(expr, Expr):
        expr = expr.expr
    else:
        expr = sympify(expr)

    # SymPy ztransform barfs on Piecewise but unilateral ZT ignores expr
    # for n < 0 so remove Piecewise.
    expr = expr.replace(var, n)        
    if expr.is_Piecewise and expr.args[0].args[1].has(n >= 0):
        expr = expr.args[0].args[0]

    expr = sym.expand(expr)        
    terms = expr.as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += ztransform_term(term, n, z)
    except ValueError:
        raise

    result = result.simplify()
    ztransform_cache[key] = result
    return result


def inverse_ztransform_ratfun(expr, z, n, **assumptions):

    damping = assumptions.get('damping', None)

    invz = symsymbol('invz', real=False)
    expr = expr.subs(z, 1 / invz)
    
    zexpr = Ratfun(expr, invz)

    Q, M, D, delay, undef = zexpr.as_QMD()

    result1 = sym.S.Zero

    if Q:
        Qpoly = sym.Poly(Q, invz)        
        C = Qpoly.all_coeffs()
        for m, c in enumerate(C):
            result1 += c * unitimpulse(n - len(C) + m + 1)

    expr = M / D
    for factor in expr.as_ordered_factors():
        if factor == sym.oo:
            return factor

    zexpr = Ratfun(expr, invz)
    poles = zexpr.poles(damping=damping)
    polesdict = {}
    for pole in poles:
        polesdict[pole.expr] = pole.n
    
    result2 = sym.S.Zero

    for pole in poles:

        p = pole.expr

        # Number of occurrences of the pole.
        o = polesdict[p]        

        if o == 0:
            continue

        if o == 1:
            r = zexpr.residue(p, poles)

            # TODO combine conjugates to get a real result
            
            result2 -= (r * p ** (-n - 1)).simplify()
            continue

        # Handle repeated poles.
        expr2 = expr * (invz - p) ** o
        for i in range(1, o + 1):
            m = o - i
            r = sym.limit(
                sym.diff(expr2, invz, m), invz, p) / sym.factorial(m)
            result2 += r * (p ** n) * n**(i - 1)

    # result1 is a sum of Dirac deltas and its derivatives so is known
    # to be causal.

    return result1, result2


def dummyvar(intnum=0):
    if intnum == 0:
        return sympify('m')
    else:
        return sympify('m_%d' % intnum)    


def inverse_ztransform_product(expr, z, n, **assumptions):

    # Handle expressions with a function of z, e.g., V(z) * Y(z), V(z)
    # / z etc.

    if assumptions.get('causal', False):
        # Assume that all functions are causal in the expression.
        n1 = sym.S.Zero
        n2 = n
    else:
        n1 = -sym.oo
        n2 = sym.oo        
    
    const, expr = factor_const(expr, z)

    factors = expr.as_ordered_factors()
    if len(factors) < 2:
        raise ValueError('Expression does not have multiple factors: %s' % expr)

    # TODO, is this useful?
    if (len(factors) > 2 and not
        isinstance(factors[1], sym.function.AppliedUndef) and
        isinstance(factors[2], sym.function.AppliedUndef)):
        factors = [factors[0], factors[2], factors[1]] + factors[3:]

    # TODO, is this useful?        
    if isinstance(factors[1], sym.function.AppliedUndef):
        terms = factors[0].as_ordered_terms()
        if len(terms) >= 2:
            result = sym.S.Zero
            for term in terms:
                result += inverse_ztransform_product(factors[1] * term, z, n)
            return result * const

    result1, result2 = inverse_ztransform_term1(factors[0], z, n)
    result = result1 + result2

    intnum = 0
    for m in range(len(factors) - 1):
        if m == 0 and isinstance(factors[1], sym.function.AppliedUndef):
            # Note, as_ordered_factors puts powers of z before the functions.
            if factors[0] == z:
                # Handle time-advance
                # Convert z * V(z) to v[n + 1]
                # TODO, fix for unilateral ZT 
                result = ztransform_func(factors[1], z, n, True)            
                result = result.subs(n, n + 1)
                continue
            elif factors[0].is_Pow and factors[0].args[0] == z and factors[0].args[1] > 0:
                # Handle higher order advances
                # Convert z ** k * V(z) to v[n + k]
                result = ztransform_func(factors[1], z, n, True)            
                result = result.subs(n, n + factors[0].args[1])
                continue                
            elif factors[0].is_Pow and factors[0].args[0] == z and factors[0].args[1] < 0:
                # Handle time-delay  1 / z ** k * V(z)
                result = ztransform_func(factors[1], z, nau, True)
                result = result.subs(n, n + factors[0].args[1])                
                continue                
        # Convert product to convolution
        dummy = dummyvar(intnum)
        intnum += 1
        result1, result2 = inverse_ztransform_term1(factors[m + 1], z, n)
        expr2 = result1 + result2
        result = sym.Sum(result.subs(n, n - dummy) * expr2.subs(n, dummy),
                         (dummy, n1, n2))
    
    return result * const


def inverse_ztransform_power(expr, z, n, **assumptions):

    # Handle expressions with a power of z.
    if not (expr.is_Pow and expr.args[0] == z):
        raise ValueError('Expression %s is not a power of z' % expr)
    exponent = expr.args[1]

    if exponent.is_positive:
        # TODO, fix
        return unitimpulse(n + exponent)

    if exponent.is_negative:
        return unitimpulse(n + exponent)

    raise ValueError('Cannot determine sign of exponent for %s' % expr)


def inverse_ztransform_sympy(expr, z, n):

    # This barfs when needing to generate Dirac deltas
    from sympy.integrals.transforms import inverse_ztransform_transform
    result = inverse_ztransform_transform(expr, n, z)
    
    if result.has(sym.InverseZtransformTransform):
        raise ValueError('Cannot determine inverse z-transform'
                         ' transform of %s with sympy' % expr)
    return result


def inverse_ztransform_term1(expr, z, n, **assumptions):

    const, expr = factor_const(expr, z)

    if expr == 1:
        return const * unitimpulse(n), sym.S.Zero
    
    if isinstance(expr, sym.function.AppliedUndef):
        # Handle V(z), 3 * V(z) etc.  If causal is True it is assumed
        # that the unknown functions are causal.  Note ztransform_func
        # just changes the name so it works as inverse_ztransform_func.
        result = ztransform_func(expr, z, n, True)
        return result * const, sym.S.Zero
    
    if expr.has(sym.function.AppliedUndef):
        return const * inverse_ztransform_product(expr, z, n,
                                               **assumptions), sym.S.Zero

    try:
        # This is the common case.
        result1, result2 = inverse_ztransform_ratfun(expr, z, n, **assumptions)
        return const * result1, const * result2
    except:
        pass

    try:
        return sym.S.Zero, const * inverse_ztransform_sympy(expr, z, n, **assumptions)
    except:
        pass

    if expr.is_Pow and expr.args[0] == z:
        return sym.S.Zero, const * inverse_ztransform_power(expr, z, n)
    
    # As last resort see if can convert to convolutions...
    return sym.S.Zero, const * inverse_ztransform_product(expr, z, n)
    
    
def inverse_ztransform_term(expr, z, n, **assumptions):

    result1, result2 = inverse_ztransform_term1(expr, z, n, **assumptions)

    if assumptions.get('causal', False):
        result2 = result2 * sym.Heaviside(n)
    
    return result1, result2


def inverse_ztransform_by_terms(expr, z, n, **assumptions):

    terms = expr.as_ordered_terms()

    result1 = sym.S.Zero
    result2 = sym.S.Zero    

    for term in terms:
        part1, part2 = inverse_ztransform_term(term, z, n, **assumptions)
        result1 += part1
        result2 += part2        
    return result1, result2


def inverse_ztransform(expr, z, n, **assumptions):
    """Calculate inverse z-transform of X(z) and return x[n].

    The unilateral z-transform cannot determine x[n] for n < 0
    unless given additional information in the way of assumptions.

    The assumptions are:
    dc -- x[n] = constant so X(z) must have the form constant * z / (z - 1)
    causal -- x[n] = 0 for n < 0.
    ac -- x[n] = A cos(a * t) + B * sin(b * t)
    """

    # TODO, simplify
    key = (expr, z, n, assumptions.get('dc', False),
           assumptions.get('ac', False),
           assumptions.get('causal', False),
           assumptions.get('damping', None))
    
    if key in inverse_ztransform_cache:
        return inverse_ztransform_cache[key]

    if expr.has(n):
        raise ValueError('Cannot inverse z-transform  %s that depends on %s' % (expr, t))
    
    if assumptions.get('dc', False):
        result = expr * (z - 1) / z
            
        free_symbols = set([symbol.name for symbol in result.free_symbols])
        if 'z' in free_symbols:
            raise ValueError('Something wonky going on, expecting dc.')
        return result

    if expr.is_Add:
        result1, result2 = inverse_ztransform_by_terms(expr, z, n, **assumptions)
    else:
        try:
            result1, result2 = inverse_ztransform_term(expr, z, n, **assumptions)
        except:
            expr = sym.expand(expr)
            result1, result2 = inverse_ztransform_by_terms(expr, z, n,
                                                           **assumptions)
        
    # result1 is known to be causal, result2 is unsure

    if assumptions.get('ac', False):
        if result1 != 0:
            raise ValueError('Inverse z-transform weirdness for %s'
                             ' with is_ac True' % expr)
        # TODO, perform more checking of the result.

    elif not assumptions.get('causal', False) and  result2 != 0:
        result2 = sym.Piecewise((result2, n >= 0))        

    result = result1 + result2    
    inverse_ztransform_cache[key] = result
    return result

from .expr import Expr

"""This module provides support for z transforms.

Copyright 2020 Michael Hayes, UCECE

"""

from .ratfun import Ratfun
from .sym import sympify, simplify, symsymbol, AppliedUndef
from .utils import factor_const, scale_shift
from .functions import UnitImpulse, UnitStep
import sympy as sym

__all__ = ('ZT', 'IZT')

ztransform_cache = {}
inverse_ztransform_cache = {}


def ztransform_func(expr, n, z, inverse=False):

    if not isinstance(expr, AppliedUndef):
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


def ztransform_sum(expr, n, z):

    const, expr = factor_const(expr, n)

    if len(expr.args) != 2:
        raise ValueError('Cannot compute z-transform of %s' % expr)

    if not isinstance(expr, sym.Sum):
        raise ValueError('Cannot compute z-transform of %s' % expr)

    # Look for integration of function
    if (isinstance(expr.args[0], AppliedUndef)
        and expr.args[1][0] == expr.args[0].args[0]
        and expr.args[1][1] == -sym.oo
        and expr.args[1][2] == n):
        return ztransform_func(expr.args[0].subs(expr.args[0].args[0], n), n, z) / (1 - 1 / z)
    
    # Look for convolution sum
    
    var = expr.args[1][0]
    if (expr.args[1][1] != -sym.oo) or (expr.args[1][2] != sym.oo):
        raise ValueError('Need indefinite limits for %s' % expr)
    
    const2, expr = factor_const(expr.args[0], n)
    if ((len(expr.args) != 2)
        or (not isinstance(expr.args[0], AppliedUndef))
        or (not isinstance(expr.args[1], AppliedUndef))):
        raise ValueError('Need sum of two functions: %s' % expr)        

    f1 = expr.args[0]
    f2 = expr.args[1]    

    # TODO: apply similarity theorem if have f(a tau) etc.
    
    if ((f1.args[0] != var or f2.args[0] != n - var)
        and (f2.args[0] != var or f1.args[0] != n - var)):
        raise ValueError('Cannot recognise convolution: %s' % expr)

    zsym = sympify(str(z))
    
    name = f1.func.__name__
    func1 = name[0].upper() + name[1:] + '(%s)' % str(zsym)

    name = f2.func.__name__
    func2 = name[0].upper() + name[1:] + '(%s)' % str(zsym)    

    F1 = sympify(func1).subs(zsym, z)
    F2 = sympify(func2).subs(zsym, z)
    
    return F1 * F2


def remove_heaviside(expr, n):

    rest = sym.S.One
    for factor in expr.as_ordered_factors():
        if (factor.is_Function and factor.func in (sym.Heaviside, UnitStep) and
            factor.args[0] == n):
            # Could remove Heaviside[n+m] where m > 0
            pass
        else:
            rest *= factor
    return rest


def ztransform_term(expr, n, z):

    const, expr = factor_const(expr, n)

    expr = remove_heaviside(expr, n)
    
    if expr.has(sym.Sum):
        try:
            return ztransform_sum(expr, n, z) * const
        except:
            pass

    nsym = sympify(str(n))
    expr = expr.replace(nsym, n)

    # foo = 1 / (1 - 1 / z)
    # factors = expr.as_ordered_factors()
    # if foo in factors:
    #     could remove factor, find ZT of rest, and then integrate....
    
    if expr.has(AppliedUndef):

        rest = sym.S.One
        expr = expr.cancel()
        for factor in expr.as_ordered_factors():
            if isinstance(factor, AppliedUndef):
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
    
    if expr.is_Function and expr.func == UnitImpulse:
        if args[0] is n:
            result = 1
        delay = n - args[0]
        if not delay.has(n):
            result = invz ** delay

    elif expr == 1:
        # Unilateral ZT
        result = 1 / (1 - invz)
        
    elif expr.is_Function and expr.func in (sym.Heaviside, UnitStep):
        if args[0] is n:
            result = 1 / (1 - invz)
        else:
            delay = n - args[0]
            if not delay.has(n):
                result = invz ** delay * 1 / (1 - invz)

    # n**i   i >= 1    
    elif (expr == n or (expr.is_Pow and args[0] == n and args[1].is_Integer and args[1]>=1) ):
        ii = 1
        # check higher order
        try:
            ii = args[1]
        except:
            pass        
        # use rule n*x[n]  o--o  -z*d/dz X(z) multiple times    
        result = 1 / (1 - invz)
        for l in range(ii):
            result = -z * sym.diff(result, z)
            result = sym.simplify(result)    
    
    # sin(b*n+c)    
    elif (expr.is_Function and expr.func == sym.sin and (args[0].as_poly(n)).is_linear):
        bb = args[0].coeff(n, 1)
        cc = args[0].coeff(n, 0)        
        result = (sym.sin(cc) + sym.sin(bb - cc) * invz) / (1 - 2 * sym.cos(bb) * invz + invz ** 2)

    # cos(b*n+c)    
    elif (expr.is_Function and expr.func == sym.cos and (args[0].as_poly(n)).is_linear):
        bb = args[0].coeff(n, 1)
        cc = args[0].coeff(n, 0)        
        result = (sym.cos(cc) - sym.cos(bb - cc) * invz) / (1 - 2 * sym.cos(bb) * invz + invz ** 2)

    # exp(b*n+c)    
    elif (expr.is_Function and expr.func == sym.exp and (args[0].as_poly(n)).is_linear):
        bb = args[0].coeff(n, 1)
        cc = args[0].coeff(n, 0)
        result = sym.exp(cc) * z / (z - sym.exp(bb))        
    
    # a**(b*n+c)
    elif (expr.is_Pow and (args[1].as_poly(n)).is_linear and expr.args[0] != n):
        bb = args[1].coeff(n, 1)
        cc = args[1].coeff(n, 0)        
        result =  args[0]**cc / (1 - args[0]**bb * invz)    

    # n * a**(b*n+c) 
    elif (expr.is_Mul and len(expr.args) == 2 and args[0] == n and args[1].is_Pow and (args[1].args[1].as_poly(n)).is_linear ): 
        # exponential part
        bb = args[1].args[1].coeff(n, 1)
        cc = args[1].args[1].coeff(n, 0)
        base_a = (args[1].args[0])    
        # result  
        result = base_a**(cc + bb) * invz / (1 - base_a**bb * invz)**2
    
       
    # n**i * a**(b*n+c)  for i>=2 
    elif (expr.is_Mul and len(expr.args) == 2 and args[0].is_Pow and (args[0].args[1].as_poly(n)).is_linear and 
           args[1].is_Pow and args[1].args[0] == n and args[1].args[1].is_Integer and args[1].args[1]>=1 ):
        ii = args[1].args[1]        
        # a**() part
        bb = args[0].args[1].coeff(n, 1)
        cc = args[0].args[1].coeff(n, 0)
        base_a = args[0].args[0]    
            
        # use rule n*x[n]  o--o  -z*d/dz X(z) multiple times    
        result = 1 / (1 - base_a**bb * invz)
        for l in range(ii):
            result = -z * sym.diff(result, z)
            result = sym.simplify(result)
        
        result *= base_a**cc    
    
    # a**(b*n+c) * sin(d*n+e)  OR  a**(b*n+c) * cos(d*n+e)
    elif (expr.is_Mul and len(expr.args) == 2 and args[0].is_Pow and ((args[0].args[1]).as_poly(n)).is_linear and args[0].args[0] != n and
           args[1].is_Function and (args[1].func == sym.sin or args[1].func == sym.cos) and 
           ((args[1].args[0]).as_poly(n)).is_linear ):
        # values for a**() part
        base_a=args[0].args[0]
        bb = args[0].args[1].coeff(n, 1)
        cc = args[0].args[1].coeff(n, 0)
        # values for sin() part
        dd = args[1].args[0].coeff(n, 1)
        ee = args[1].args[0].coeff(n, 0)
        
        pre_sign = 1
        if (expr.args[1]).func == sym.cos:
            pre_sign = -1  
        
        result = base_a**cc * (expr.args[1].func(ee) + pre_sign * base_a**bb * expr.args[1].func(dd - ee) * invz) / (1 - 2 * sym.cos(dd) * base_a**bb * invz + (base_a**bb * invz) ** 2)
    
    
    # n * sin(b*n+c)  OR  n * cos(b*n+c) for i>=1
    elif (expr.is_Mul and len(expr.args) == 2 and args[0] == n and 
          args[1].is_Function and (args[1].func == sym.sin or args[1].func == sym.cos) 
          and ((expr.args[1].args[0]).as_poly(n)).is_linear ) :
        # find parameter of sin() part
        bb = (args[1].args[0]).coeff(n, 1)
        cc = (args[1].args[0]).coeff(n, 0)        

        pre_sign = 1
        if args[1].func == sym.cos:
            pre_sign = -1        

        result = (args[1].func(cc) + pre_sign * args[1].func(bb - cc) * invz) / (1 - 2 * sym.cos(bb) * invz + invz ** 2)    
        result = -z * sym.diff(result, z)
            
    
    
        
    # n**i * sin(b*n+c)  OR  n**i * cos(b*n+c) for i>=1
    elif (expr.is_Mul and len(expr.args) == 2 and args[0].is_Pow and args[0].args[0] == n and args[0].args[1].is_Integer and args[0].args[1]>=1 and 
          args[1].is_Function and (args[1].func == sym.sin or args[1].func == sym.cos) 
             and ((expr.args[1].args[0]).as_poly(n)).is_linear ) :
        ii = args[0].args[1]
        # find parameter of sin() part
        bb = (args[1].args[0]).coeff(n, 1)
        cc = (args[1].args[0]).coeff(n, 0)        
        
        pre_sign = 1
        if args[1].func == sym.cos:
            pre_sign = -1        
        
        result = (args[1].func(cc) + pre_sign * args[1].func(bb - cc) * invz) / (1 - 2 * sym.cos(bb) * invz + invz ** 2)    
        for l in range(ii):
            result = -z * sym.diff(result, z)
            result = sym.simplify(result)  
    
    # n * a**(b*n+c) * sin(d*n+e)  OR n * a**(b*n+c) * cos(d*n+e) 
    elif (expr.is_Mul and len(expr.args) == 3 and args[0] == n and args[2].is_Function and (args[2].func == sym.sin or args[2].func == sym.cos) 
           and ((args[2].args[0]).as_poly(n)).is_linear and
           args[1].is_Pow and ((args[1].args[1]).as_poly(n)).is_linear and args[1].args[0] != n ) :
        # find parameter of sin(), cos() part
        dd = (args[2].args[0]).coeff(n, 1)
        ee = (args[2].args[0]).coeff(n, 0) 
        # values for a**() part
        base_a = args[1].args[0]
        bb = args[1].args[1].coeff(n, 1)
        cc = args[1].args[1].coeff(n, 0)        
        # check for sign in numerator
        pre_sign = 1
        if args[2].func == sym.cos:
            pre_sign = -1        
        # result
        result = (args[2].func(ee) + pre_sign * args[2].func(dd - ee) * base_a**bb * invz) / (1 - 2 * sym.cos(dd) * base_a**bb * invz + (base_a**bb*invz) ** 2)         
        result = -z*sym.diff(result, z)
        result = sym.simplify(result)
    
    
    # n**i * a**(b*n+c) * sin(d*n+e)  OR  n**i * a**(b*n+c) * cos(d*n+e)  i>=2
    elif (expr.is_Mul and len(expr.args) == 3 and args[1].is_Pow and args[1].args[0] == n and args[1].args[1].is_Integer and args[1].args[1]>=2 and 
          args[2].is_Function and (args[2].func == sym.sin or args[2].func == sym.cos) and ((args[2].args[0]).as_poly(n)).is_linear and
           args[0].is_Pow and ((args[0].args[1]).as_poly(n)).is_linear and args[0].args[0] != n ):  
        ii = args[1].args[1]
        # find parameter of sin(), cos() part
        dd = (args[2].args[0]).coeff(n, 1)
        ee = (args[2].args[0]).coeff(n, 0) 
        # values for a**() part
        base_a = expr.args[0].args[0]
        bb = args[0].args[1].coeff(n, 1)
        cc = args[0].args[1].coeff(n, 0)  
        # check for sign in numerator
        pre_sign = 1
        if args[2].func == sym.cos:
            pre_sign = -1
            
        # use rule n*x[n]  o--o  -z*d/dz X(z) multiple times    
        result = (args[2].func(ee) + pre_sign * args[2].func(dd - ee) * base_a**bb * invz) / (1 - 2 * sym.cos(dd) * base_a**bb * invz + (base_a**bb * invz) ** 2)    
        for l in range(ii):
            result = -z * sym.diff(result, z)
            result = sym.simplify(result)        
        result *= base_a**cc


    # n * exp(a * n)
    elif (expr.is_Mul and len(args) == 2 and args[0] == n and
          args[1].is_Function and args[1].func == sym.exp):
        aconst, aexpr = factor_const(args[1].args[0], n)
        if aexpr == n:
            result = z * sym.exp(aconst) / (z - sym.exp(aconst))**2            

    
    # n**i * exp(a * n)   ;   i>=2
    elif (expr.is_Mul and len(expr.args) == 2 and args[1].is_Function and args[1].func == sym.exp and 
          args[0].is_Pow and args[0].args[0] == n and args[0].args[1].is_Integer and args[0].args[1]>=1 ):
        ii = args[0].args[1]        
        # exp**() part
        aconst, aexpr = factor_const(args[1].args[0], n)    

        # use rule n*x[n]  o--o  -z*d/dz X(z) multiple times    
        result = 1 / (1 - sym.exp(aconst) * invz)
        for l in range(ii):
            result = -z * sym.diff(result, z)
            result = sym.simplify(result)   


    # sin(a * n + b) * exp(c * n + d) and exp(c * n + d) * sin(a * n + b)
    elif (expr.is_Mul and args[0].is_Function and args[1].is_Function and
          args[0].func == sym.exp and args[1].func == sym.sin and ((args[1].args[0]).as_poly(n)).is_linear and ((args[0].args[0]).as_poly(n)).is_linear ):
        # find parameter of sin() part
        aa = (args[1].args[0]).coeff(n, 1)
        bb = (args[1].args[0]).coeff(n, 0)
        # exp part
        cc = (args[0].args[0]).coeff(n, 1)
        dd = (args[0].args[0]).coeff(n, 0)        
        # result
        num = sym.sin(bb) + sym.exp(cc) * sym.sin(aa-bb) * invz
        den = (1 - 2 * invz * sym.exp(cc) * sym.cos(aa) + sym.exp(2 * cc)*invz**2)
        result = sym.exp(dd) * num / den


    # cos(a * n + b) * exp(c * n + d) and exp(c * n + d) * cos(a * n + b)
    elif (expr.is_Mul and args[0].is_Function and args[1].is_Function and
          args[1].func == sym.exp and args[0].func == sym.cos and ((args[1].args[0]).as_poly(n)).is_linear and ((args[0].args[0]).as_poly(n)).is_linear ):
        # find parameter of cos() part
        aa = (args[0].args[0]).coeff(n, 1)
        bb = (args[0].args[0]).coeff(n, 0)
        # exp part
        cc = (args[1].args[0]).coeff(n, 1)
        dd = (args[1].args[0]).coeff(n, 0)        
        # result
        num = sym.cos(bb) - sym.exp(cc) * sym.cos(aa-bb) * invz
        den = (1 - 2 * invz * sym.exp(cc) * sym.cos(aa) + sym.exp(2 * cc)*invz**2)
        result = sym.exp(dd) * num / den            
    

    
    # n * cos(a*n+b) * exp(c*n+d) 
    elif (expr.is_Mul and len(expr.args) == 3 and args[0] == n and args[1].is_Function and args[1].func == sym.cos and ((args[1].args[0]).as_poly(n)).is_linear and
           args[2].is_Function and args[2].func == sym.exp and ((args[2].args[0]).as_poly(n)).is_linear ) :
        # find parameter of cos() part
        aa = (args[1].args[0]).coeff(n, 1)
        bb = (args[1].args[0]).coeff(n, 0) 
        # values for exp() part
        cc = args[2].args[0].coeff(n, 1)
        dd = args[2].args[0].coeff(n, 0)        
        # result
        result = (sym.cos(bb) - sym.cos(aa-bb) * sym.exp(cc) * invz) / (1 - 2 * sym.cos(aa) * sym.exp(cc) * invz + (sym.exp(cc)*invz) ** 2)         
        result = -z*sym.diff(result, z)
        result = sym.simplify(result)    
        result*=sym.exp(dd)
        
    
    # n**ii *  cos(a*n+b) * exp(c*n+d) 
    elif (expr.is_Mul and len(expr.args) == 3 and args[0].is_Pow and args[0].args[0] == n and args[0].args[1].is_Integer and args[0].args[1]>=1 and 
          args[1].is_Function and args[1].func == sym.cos and ((args[1].args[0]).as_poly(n)).is_linear and
          args[2].is_Function and args[2].func == sym.exp and ((args[2].args[0]).as_poly(n)).is_linear ) :
        # find parameter of cos() part
        aa = (args[1].args[0]).coeff(n, 1)
        bb = (args[1].args[0]).coeff(n, 0) 
        # values for exp() part
        cc = args[2].args[0].coeff(n, 1)
        dd = args[2].args[0].coeff(n, 0)  
        # ii
        ii = args[0].args[1]
        # result
        result = (sym.cos(bb) - sym.cos(aa-bb) * sym.exp(cc) * invz) / (1 - 2 * sym.cos(aa) * sym.exp(cc) * invz + (sym.exp(cc)*invz) ** 2)         
        for l in range(ii):
            result = -z * sym.diff(result, z)
            result = sym.simplify(result)   
        result*=sym.exp(dd)    
    
    # n  * exp(c*n+d) *  sin(a*n+b)
    elif (expr.is_Mul and len(expr.args) == 3 and args[0] == n and args[2].is_Function and args[2].func == sym.sin and ((args[2].args[0]).as_poly(n)).is_linear and
          args[1].is_Function and args[1].func == sym.exp and ((args[1].args[0]).as_poly(n)).is_linear ) :
        # find parameter of sin() part
        aa = (args[2].args[0]).coeff(n, 1)
        bb = (args[2].args[0]).coeff(n, 0) 
        # values for exp() part
        cc = args[1].args[0].coeff(n, 1)
        dd = args[1].args[0].coeff(n, 0)        
        # result
        result = (sym.sin(bb) + sym.sin(aa-bb) * sym.exp(cc) * invz) / (1 - 2 * sym.cos(aa) * sym.exp(cc) * invz + sym.exp(2*cc)* invz ** 2)         
        result = -z*sym.diff(result, z)
        result = sym.simplify(result) 
        result*=sym.exp(dd)
    
    
    # n**ii *  sin(a*n+b) * exp(c*n+d) 
    elif (expr.is_Mul and len(expr.args) == 3 and args[0].is_Pow and args[0].args[0] == n and args[0].args[1].is_Integer and args[0].args[1]>=1 and 
          args[2].is_Function and args[2].func == sym.sin and ((args[2].args[0]).as_poly(n)).is_linear and
          args[1].is_Function and args[1].func == sym.exp and ((args[1].args[0]).as_poly(n)).is_linear ) :
        # find parameter of sin() part
        aa = (args[2].args[0]).coeff(n, 1)
        bb = (args[2].args[0]).coeff(n, 0) 
        # values for exp() part
        cc = args[1].args[0].coeff(n, 1)
        dd = args[1].args[0].coeff(n, 0)  
        # ii
        ii = args[0].args[1]
        # result
        result = (sym.sin(bb) + sym.sin(aa-bb) * sym.exp(cc) * invz) / (1 - 2 * sym.cos(aa) * sym.exp(cc) * invz + (sym.exp(cc)*invz) ** 2)         
        for l in range(ii):
            result = -z * sym.diff(result, z)
            result = sym.simplify(result)    
        result*=sym.exp(dd)
    
    
    # exp(a * n) * exp(b * n)
    elif (expr.is_Mul and args[0].is_Function and args[1].is_Function and
          args[0].func == sym.exp and args[1].func == sym.exp):
        aconst, aexpr = factor_const(args[1].args[0], n)
        bconst, bexpr = factor_const(args[0].args[0], n)
        if aexpr == n and bexpr == n:
            result = z / (z - sym.exp(aconst + bconst))    
    
        
    if result is None:
        # Use m instead of n to avoid n and z in same expr.
        # TODO, check if m already used...
        msym = sympify('m', real=True)        
        nsym = sympify(str(n))        
        zsym = sympify(str(z))
        result = sym.Sum(expr.subs(nsym, msym) * zsym**msym, (msym, 0, sym.oo))
        
    return const * result


def ztransform(expr, n, z, evaluate=True):
    """Compute unilateral z-transform of expr with lower index 0.

    Undefined functions such as v[n] are converted to V(z)

    """

    if expr.is_Equality:
        return sym.Eq(ztransform(expr.args[0], n, z),
                      ztransform(expr.args[1], n, z))

    if not evaluate:
        result = sym.Sum(expr * z**(-n), (n, 0, sym.oo))
        return result
    
    const, expr = factor_const(expr, n)    
    key = (expr, n, z)
    if key in ztransform_cache:
        return const * ztransform_cache[key]

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
    return (const * result).simplify()


def inverse_ztransform_ratfun(expr, z, n, **assumptions):

    expr = expr / z
    
    # Handle special case 1 / (z**m * (z - 1)) since this becomes u[n - m]
    # The default method produces u[n] - delta[n] for u[n-1].  This is correct
    # but can be simplified.
    # In general, 1 / (z**m * (z - a)) becomes a**n * u[n - m]

    if (len(expr.args) == 2 and expr.args[1].is_Pow and
        expr.args[1].args[0].is_Add and
        expr.args[1].args[0].args[0] == -1 and
        expr.args[1].args[0].args[1] == z):

        delay = None
        if expr.args[0] == z:
            delay = 1
        elif expr.args[0].is_Pow and expr.args[0].args[0] == z:
            a = expr.args[0].args[1]
            if a.is_positive:
                print('Warning, dodgy z-transform.  Have advance of unit step.')
            elif not a.is_negative:
                print('Warning, dodgy z-transform.  May have advance of unit step.')                
            delay = -a
        elif (expr.args[0].is_Pow and expr.args[0].args[0].is_Pow and
              expr.args[0].args[0].args[0] == z and
              expr.args[0].args[0].args[1] == -1):              
            a = expr.args[0].args[1]
            if a.is_negative:
                print('Warning, dodgy z-transform.  Have advance of unit step.')
            elif not a.is_positive:
                print('Warning, dodgy z-transform.  May have advance of unit step.')                
            delay = a            

        if delay is not None:
            return UnitStep(n - delay), sym.S.Zero
        
    damping = assumptions.get('damping', None)
    
    zexpr = Ratfun(expr, z)

    Q, M, D, delay, undef = zexpr.as_QMD()

    cresult = sym.S.Zero
    uresult = sym.S.Zero

    if Q:
        Qpoly = sym.Poly(Q, z)        
        C = Qpoly.all_coeffs()
        for m, c in enumerate(C):
            cresult += c * UnitImpulse(n - len(C) + m + 1)

    # There is problem with determining residues if
    # have 1/(z*(-a/z + 1)) instead of 1/(-a + z).  Hopefully,
    # simplify will fix things...
    expr = (M / D).simplify()
    for factor in expr.as_ordered_factors():
        if factor == sym.oo:
            return factor, factor

    zexpr = Ratfun(expr, z)
    poles = zexpr.poles(damping=damping)
    polesdict = {}
    for pole in poles:
        polesdict[pole.expr] = pole.n
    
    for pole in poles:

        p = pole.expr

        # Number of occurrences of the pole.
        o = polesdict[p]        

        if o == 0:
            continue

        if o == 1:
            r = zexpr.residue(p, poles)

            # TODO combine conjugates to get a real result
            # See laplace.py

            if p == 0:
                cresult += r * UnitImpulse(n)
            else:
                uresult += r * p ** n
            continue

        # Handle repeated poles.
        
        expr2 = expr * (z - p) ** o
        expr2 = expr2.simplify()
        bino = 1
        for i in range(1, o + 1):
            m = o - i
            r = sym.limit(
                sym.diff(expr2, z, m), z, p) / sym.factorial(m)

            if p == 0:
                cresult += r * UnitImpulse(n - i + 1)
            else:            
                uresult += r * bino * p **(n - i + 1) / sym.factorial(i - 1)
                bino *= n - i + 1

    # cresult is a sum of Dirac deltas and its derivatives so is known
    # to be causal.

    return cresult, uresult


def dummyvar(intnum=0):
    if intnum == 0:
        return sympify('m', real=True)
    else:
        return sympify('m_%d' % intnum, real=True)


def inverse_ztransform_product(expr, z, n, **assumptions):

    # Handle expressions with a function of z, e.g., V(z) * Y(z), V(z)
    # / z etc.

    zsym = sympify(str(z))    

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

    if (len(factors) == 3 and factors[0] == z and
        isinstance(factors[2], AppliedUndef) and
        factors[1].is_Pow and factors[1].args[1] == -1 and
        factors[1].args[0].is_Add and factors[1].args[0].args[0] == -1
        and factors[1].args[0].args[1] == z):
        # Handle cumulative sum  z / (z - 1) * V(z)
        m = dummyvar()
        result = ztransform_func(factors[2], z, m, True)                
        result = sym.Sum(result, (m, n1, n))
        return result
    
    # TODO, is this useful?
    if (len(factors) > 2 and not
        isinstance(factors[1], AppliedUndef) and
        isinstance(factors[2], AppliedUndef)):
        factors = [factors[0], factors[2], factors[1]] + factors[3:]

    # TODO, is this useful?        
    if isinstance(factors[1], AppliedUndef):
        terms = factors[0].as_ordered_terms()
        if len(terms) >= 2:
            result = sym.S.Zero
            for term in terms:
                result += inverse_ztransform_product(factors[1] * term, z, n)
            return const * result

    cresult, uresult = inverse_ztransform_term1(factors[0], z, n)
    result = cresult + uresult

    intnum = 0
    for m in range(len(factors) - 1):
        if m == 0 and isinstance(factors[1], AppliedUndef):
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
                result = ztransform_func(factors[1], z, n, True)
                result = result.subs(n, n + factors[0].args[1])                
                continue
            elif (factors[0].is_Pow and
                  factors[0].args[0].is_Add and
                  factors[0].args[1] == -1 and                  
                  factors[0].args[0].args[0] == 1 and
                  factors[0].args[0].args[1].is_Mul and
                  factors[0].args[0].args[1].args[0] == -1 and
                  factors[0].args[0].args[1].args[1].is_Pow and
                  factors[0].args[0].args[1].args[1].args[1] == -1 and
                  factors[0].args[0].args[1].args[1].args[0] is zsym):
                # Handle cumulative sum  1 / (1 - 1 / z) * V(z)
                m = dummyvar(intnum)
                result = ztransform_func(factors[1], z, m, True)                
                intnum += 1
                result = sym.Sum(result, (m, n1, n))
                continue

        # Convert product to convolution
        dummy = dummyvar(intnum)
        intnum += 1
        cresult, uresult = inverse_ztransform_term1(factors[m + 1], z, n)
        expr2 = cresult + uresult
        kernel = result.subs(n, n - dummy) * expr2.subs(n, dummy)
        sresult = sym.Sum(kernel, (dummy, n1, n2))
        result = sresult
    
    return const * result


def inverse_ztransform_power(expr, z, n, **assumptions):

    # Handle expressions with a power of z.
    if (expr.is_Pow and expr.args[0] == z):
        exponent = expr.args[1]

        if exponent.is_positive:
            print('Warning, dodgy z-transform.  Have advance of unit impulse.')
        elif not exponent.is_negative:
            print('Warning, dodgy z-transform.  May have advance of unit impulse.')
        
        return UnitImpulse(n + exponent), sym.S.Zero

    # Handle expressions with a power of (1 / z).
    if (expr.is_Pow and expr.args[0].is_Pow and
        expr.args[0].args[0] == z and expr.args[0].args[1] == -1):

        exponent = expr.args[1]

        if exponent.is_negative:
            print('Warning, dodgy z-transform.  Have advance of unit impulse.')
        elif not exponent.is_positive:
            print('Warning, dodgy z-transform.  May have advance of unit impulse.')
        
        return UnitImpulse(n - exponent), sym.S.Zero        

    raise ValueError('Expression %s is not a power of z' % expr)


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
        return const * UnitImpulse(n), sym.S.Zero
    
    if isinstance(expr, AppliedUndef):
        # Handle V(z), 3 * V(z) etc.  If causal is True it is assumed
        # that the unknown functions are causal.  Note ztransform_func
        # just changes the name so it works as inverse_ztransform_func.
        result = ztransform_func(expr, z, n, True)
        return const * result, sym.S.Zero
    
    if expr.has(AppliedUndef):
        return const * inverse_ztransform_product(expr, z, n,
                                                  **assumptions), sym.S.Zero

    if expr == z:
        print('Warning, dodgy z-transform.  Have advance of unit impulse.') 
        return const * UnitImpulse(n + 1), sym.S.Zero        

    if (expr.is_Pow and
        (expr.args[0] == z or
         (expr.args[0].is_Pow and
          expr.args[0].args[0] == z and expr.args[0].args[1] == -1))):
        cresult, uresult = inverse_ztransform_power(expr, z, n, **assumptions)
        return const * cresult, const * uresult        
        
    try:
        # This is the common case.
        cresult, uresult = inverse_ztransform_ratfun(expr, z, n, **assumptions)
        return const * cresult, const * uresult
    except:
        pass

    try:
        return sym.S.Zero, const * inverse_ztransform_sympy(expr, z, n, **assumptions)
    except:
        pass

    # As last resort see if can convert to convolutions...
    return sym.S.Zero, const * inverse_ztransform_product(expr, z, n)
    
    
def inverse_ztransform_term(expr, z, n, **assumptions):

    cresult, uresult = inverse_ztransform_term1(expr, z, n, **assumptions)

    if assumptions.get('causal', False):
        uresult = uresult * UnitStep(n)
    
    return cresult, uresult


def inverse_ztransform_by_terms(expr, z, n, **assumptions):

    terms = expr.as_ordered_terms()

    cresult = sym.S.Zero
    uresult = sym.S.Zero    

    for term in terms:
        part1, part2 = inverse_ztransform_term(term, z, n, **assumptions)
        cresult += part1
        uresult += part2        
    return cresult, uresult


def inverse_ztransform_make(n, const, cresult, uresult, **assumptions):

    result = const * (cresult + uresult)
    result = result.simplify()
    
    if assumptions.get('dc', False):
        free_symbols = set([symbol.name for symbol in result.free_symbols])
        if 'n' in free_symbols:
            raise ValueError('Something wonky going on, expecting dc.')
    
    elif assumptions.get('ac', False):

        if cresult != 0:
            raise ValueError('Inverse z-transform weirdness for %s with is_ac True' % result)
        # TODO, perform more checking of the result.
        
    elif not assumptions.get('causal', False):

        # Cannot determine result for n < 0
        result = sym.Piecewise((result, n >= 0))
            
    return result    


def inverse_ztransform1(expr, z, n, **assumptions):

    const, expr = factor_const(expr, z)
    
    key = (expr, z, n, 
           assumptions.get('causal', False),
           assumptions.get('damping', None))
    
    if key in inverse_ztransform_cache:
        cresult, uresult = inverse_ztransform_cache[key]        
        return const, cresult, uresult        

    try:
        if expr.is_Add:
            cresult, uresult = inverse_ztransform_by_terms(expr, z, n, **assumptions)
        else:
            try:
                cresult, uresult = inverse_ztransform_term(expr, z, n, **assumptions)
            except:
                expr = sym.expand(expr)
                cresult, uresult = inverse_ztransform_by_terms(expr, z, n,
                                                               **assumptions)
    except:
        raise ValueError('Cannot determine z-transform of %s' % expr)
        
    # cresult is known to be causal, uresult is unsure
    inverse_ztransform_cache[key] = cresult, uresult
    return const, cresult, uresult    


def inverse_ztransform(expr, z, n, **assumptions):
    """Calculate inverse z-transform of X(z) and return x[n].

    The unilateral z-transform cannot determine x[n] for n < 0
    unless given additional information in the way of assumptions.

    The assumptions are:
    dc -- x[n] = constant so X(z) must have the form constant * z / (z - 1)
    causal -- x[n] = 0 for n < 0.
    ac -- x[n] = A cos(a * t) + B * sin(b * t)
    """

    if expr.is_Equality:
        return sym.Eq(inverse_ztransform(expr.args[0], z, n, **assumptions),
                      inverse_ztransform(expr.args[1], z, n, **assumptions))

    if expr.has(n):
        raise ValueError('Cannot inverse z-transform  %s that depends on %s' % (expr, n))

    const, cresult, uresult = inverse_ztransform1(expr, z, n, **assumptions)
    return inverse_ztransform_make(n, const, cresult, uresult, **assumptions)


def ZT(expr, n, z, **assumptions):
    """Compute unilateral Z-Transform transform of expr with lower limit 0.

    Undefined functions such as v[n] are converted to V(z)."""
    
    return ztransform(expr, n, z, **assumptions)


def IZT(expr, z, n, **assumptions):
    """Calculate inverse z-Transform of X(s) and return x[n].

    The unilateral Z-Transform transform cannot determine x[n] for n < 0
    unless given additional information in the way of assumptions.

    The assumptions are:
    dc -- x[n] = constant
    causal -- x[n] = 0 for n < 0.
    ac -- x[n] = A cos(a * n) + B * sin(b * n)
    """
    
    return inverse_ztransform(expr, z, n, **assumptions)

from .expr import Expr


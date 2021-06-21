"""This module provides support for z transforms.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from .transformer import UnilateralForwardTransformer
from .ratfun import Ratfun
from .sym import sympify, simplify, symsymbol, AppliedUndef
from .utils import factor_const, scale_shift
from .extrafunctions import UnitImpulse, UnitStep
import sympy as sym
from sympy.simplify.fu import TR6, TR9

__all__ = ('ZT', )


# Check the structure of the given expression
def is_multiplied_with(expr, n, cmp, ret):

    ret_flag = False   
    # Check for multiplication  with n
    if cmp == 'n' and expr == n:  #only n
        ret += [n]
        ret_flag = True  
    elif (cmp == 'n' and expr.is_Pow and expr.args[0] == n and
          expr.args[1].is_integer and expr.args[1] >= 1):     # only n**i
        ret += [n]
        ret_flag = True        
    elif cmp == 'n' and expr.is_Mul:   # multiplication with n
        for i in range(len(expr.args)):
            if (expr.args[i].is_Pow and expr.args[i].args[0] == n and
                expr.args[i].args[1].is_integer and expr.args[i].args[1] >= 1):
                ret += [n]
                ret_flag = True
                break
            elif expr.args[i] == n:
                ret += [n]
                ret_flag = True
                break

    # Check for multiplication with a**(b*n+c)            
    elif (cmp == 'a**n' and expr.is_Pow and
          ((expr.args[1]).as_poly(n)).is_linear and (not expr.args[0].has(n))):
        ret += [expr]
        ret_flag = True
    elif cmp == 'a**n' and expr.is_Mul:   
        for i in range(len(expr.args)):
            if ((expr.args[i].is_Pow and
                 ((expr.args[i].args[1]).as_poly(n)).is_linear and
                 (not expr.args[i].args[0].has(n)))):
                ret += [expr.args[i]]   
                ret_flag = True
                break

    # Check for multiplication with exp(b*n+c)        
    elif (cmp == 'exp(n)' and len(expr.args) == 1 and expr.is_Function and
          expr.func == sym.exp):
        ret += [expr]
        ret_flag = True
    elif cmp == 'exp(n)' and expr.is_Mul:   
        for i in range(len(expr.args)):
            if (expr.args[i].is_Function and expr.args[i].func == sym.exp and
                ((expr.args[i].args[0]).as_poly(n)).is_linear):
                ret += [expr.args[i]]   
                ret_flag = True
                break                     

    # Check for Multiplication with u(n-n0)
    elif cmp == 'UnitStep' and expr.is_Mul:
        for i in range(len(expr.args)):
            if expr.args[i].is_Function and expr.args[i].func in (sym.Heaviside, UnitStep):
                ret += [expr.args[i]]   
                ret_flag = True
                break 
   
   # check for Multiplication with sin
    elif (cmp == 'sin(n)' and len(expr.args) == 1 and expr.is_Function and  # sin only
          expr.func == sym.sin):
        ret += [expr]
        ret_flag = True
    elif cmp == 'sin(n)' and expr.is_Mul:   
        for i in range(len(expr.args)):
            if (expr.args[i].is_Function and expr.args[i].func == sym.sin and
                ((expr.args[i].args[0]).as_poly(n)).is_linear):
                ret += [expr.args[i]]   
                ret_flag = True
                break 
            
    # check for Multiplication with cos
    elif (cmp == 'cos(n)' and len(expr.args) == 1 and expr.is_Function and  # sin only
          expr.func == sym.cos):
        ret += [expr]
        ret_flag = True
    elif cmp == 'cos(n)' and expr.is_Mul:   
        for i in range(len(expr.args)):
            if (expr.args[i].is_Function and expr.args[i].func == sym.cos and
                ((expr.args[i].args[0]).as_poly(n)).is_linear):
                ret += [expr.args[i]]   
                ret_flag = True
                break             


    return ret_flag


class ZTransformer(UnilateralForwardTransformer):

    name = 'z-transform'
    
    def noevaluate(self, expr, n, z):

        foo = expr * z**(-n)
        result = sym.Sum(foo, (n, 0, sym.oo))
        return result

    def rewrite(self, expr, var):
        # This is needed to handle expressions like (2*n + 3)**2
        return sym.expand(expr)
    
    def check(self, expr, n, z, **assumptions):

        if expr.has(z):
            self.error('Expression depends on z')

    def key(self, expr, n, z, **assumptions):
        return expr, n, z,

    def func(self, expr, n, z):

        if not isinstance(expr, AppliedUndef):
            raise ValueError('Expecting function for %s' % expr)

        scale, shift = scale_shift(expr.args[0], n)    

        zsym = sympify(str(z))

        # Convert v(n) to V(z), etc.
        name = expr.func.__name__
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


    def sum(self, expr, n, z):

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
            return self.func(expr.args[0].subs(expr.args[0].args[0], n), n, z) / (1 - 1 / z)

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

    def term(self, expr, n, z):

        const, expr = factor_const(expr, n)

        if expr.has(sym.Sum):
            try:
                return self.sum(expr, n, z) * const
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
                    result = self.func(factor, n, z)
                else:
                    if factor.has(n):
                        raise ValueError('TODO: need derivative of undefined'
                                         ' function for %s' % factor)
                    rest *= factor
            return result * rest * const

        invz = z ** -1

        result = None
        args = expr.args
        xn_fac = []

        if expr.is_Function and expr.func == UnitImpulse:
            if args[0] is n:
                result = 1
            delay = n - args[0]
            if not delay.has(n):
                result = invz ** delay

        elif expr == 1:
            # Unilateral ZT
            result = 1 / (1 - invz)

        elif (expr.is_Function and
              expr.func in (sym.Heaviside, UnitStep, sym.sign)):
            if args[0] is n:
                result = 1 / (1 - invz)
            else:
                delay = n - args[0]
                if not delay.has(n):
                    result = invz ** delay * 1 / (1 - invz)

        # sin(b*n+c)    
        elif (expr.is_Function and expr.func == sym.sin and (args[0].as_poly(n)).is_linear):
            bb = args[0].coeff(n, 1)
            cc = args[0].coeff(n, 0)        
            result =  (sym.sin(cc) + sym.sin(bb - cc) * invz) / (1 - 2 * sym.cos(bb) * invz + invz**2) 
            result = sym.simplify(result)

        # cos(b*n+c)    
        elif (expr.is_Function and expr.func == sym.cos and (args[0].as_poly(n)).is_linear):
            bb = args[0].coeff(n, 1)
            cc = args[0].coeff(n, 0)        
            result = (sym.cos(cc) - sym.cos(bb - cc) * invz) / (1 - 2 * sym.cos(bb) * invz + invz**2)  
            result = sym.simplify(result)

        # Multiplication with n       use n*x(n)  o--o  -z d/dz X(z)
        elif is_multiplied_with(expr, n, 'n', xn_fac):
            expr = expr / xn_fac[0]
            X = self.term(expr, n, z)
            result = sym.simplify(-z * sym.diff(X, z))

        # Multiplication with a**(b*n+c)   use    lam**n *x(n)  o--o  X(z/lam)
        elif is_multiplied_with(expr, n, 'a**n', xn_fac):
            expr /= xn_fac[0]
            ref = xn_fac[0].args
            lam = ref[0]
            bb = ref[1].coeff(n, 1)
            cc = ref[1].coeff(n, 0)              
            X = self.term(expr, n, z)
            result = lam**cc * sym.simplify(X.subs(z, z / lam**bb)) 

        # Multiplication with exp(b*n+c)   use    exp**n *x(n)  o--o  X(z/exp(1))
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac):
            expr /= xn_fac[0]
            ref = xn_fac[0].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)              
            X = self.term(expr, n, z)
            result = sym.exp(cc) * sym.simplify(X.subs(z, z / sym.exp(bb))) 

        # Multiplication with u(n-n0) * x(n)    o--o     X(z)- (x[0] - x[1]/z - .. - x[n0-1]/z**(n0-1)   
        elif is_multiplied_with(expr, n, 'UnitStep', xn_fac):
            expr /= xn_fac[0]
            delay = n - xn_fac[0].args[0]
            X = self.term(expr, n, z)
            sum_X = 0
            for ii in range(delay):
                sum_X -= expr.subs(n, ii) * invz**ii
            result = X + sum_X               

        if result is None:
            # Use m instead of n to avoid n and z in same expr.
            # TODO, check if m already used...
            msym = sympify('m', real=True)        
            nsym = sympify(str(n))        
            zsym = sympify(str(z))
            result = sym.Sum(expr.subs(nsym, msym) * zsym**msym, (msym, 0, sym.oo))

        return const * result


ztransformer = ZTransformer()


def ZT(expr, n, z, evaluate=True, **assumptions):
    """Compute unilateral Z-Transform transform of expr with lower limit 0.

    Undefined functions such as v[n] are converted to V(z)."""
    
    return ztransformer.transform(expr, n, z,
                                  evaluate=evaluate, **assumptions)

def ztransform(expr, n, z, evaluate=True, **assumptions):
    """Compute unilateral Z-Transform transform of expr with lower limit 0.

    Undefined functions such as v[n] are converted to V(z)."""

    return ztransformer.transform(expr, n, z,
                                  evaluate=evaluate, **assumptions)    


from .expr import Expr

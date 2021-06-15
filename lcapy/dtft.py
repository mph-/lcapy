"""This module provides support for discrete-time Fourier transforms. 

 It calculates the discrete-time Fourier transform using:

   X(f) = \sum_{n=\infty}^{\infty} x(n) e^{-j * 2 * \pi * n * dt * f}

Copyright 2021 Michael Hayes, UCECE

"""

import sympy as sym
from sympy import oo, DiracDelta
from .transformer import BilateralForwardTransformer
from .sym import sympify, AppliedUndef, j, pi
from .dsym import dt
from .utils import factor_const, scale_shift
from .sym import symsymbol
from .ztransform import is_multiplied_with
from .extrafunctions import UnitImpulse, UnitStep, sincu, sincn, rect


__all__ = ('DTFT', )


class DTFTTransformer(BilateralForwardTransformer):

    name = 'DTFT'
    
    def key(self, expr, n, f, **assumptions):
        return expr, n, f, assumptions.get('images', 0)

    def noevaluate(self, expr, n, f):

        foo = expr * sym.exp(-2 * j * pi * n * dt * f)
        result = sym.Sum(foo, (n, -oo, oo))
        return result

    def check(self, expr, n, f, images=0, **assumptions):

        self.images = images
        if images == oo:
            self.m1 = -oo
            self.m2 = oo
        else:
            self.m1 = -(images // 2)
            self.m2 = self.m1 + self.images
        
        if expr.has(f):
            self.error('Expression depends on f')
        
        if expr.is_Piecewise and expr.args[0].args[1].has(n >= 0):
            self.error('Expression is unknown for n < 0' % expr)

    def add_images(self, expr, f):
        if self.m1 == self.m2:
            return expr

        msym = symsymbol('m', integer=True)
        foo = expr.replace(f, f - msym / dt)
        result = sym.Sum(foo, (msym, self.m1, self.m2))
        return result
            
    def sympy(self, expr, n, f):

        foo = expr * sym.exp(-2 * j * pi * n * dt * f)
        result = sym.summation(foo, (n, -oo, oo))

        return result

    def func(self, expr, n, f):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], n)

        fsym = sympify(str(f))

        # Convert v(n) to V(f), etc.
        name = expr.func.__name__
        func = name[0].upper() + name[1:] + '(%s)' % f

        result = sympify(func).subs(fsym, f / scale) / abs(scale)

        if shift != 0:
            result = result * sym.exp(2 * sym.I * sym.pi * f * shift / scale)
        # Perhaps return X_(1/dt)(f) but how to denote?
        return self.add_images(result, f)

    def function(self, expr, n, f):

        # Handle expressions with a function of FOO, e.g.,
        # v(n), v(n) * y(n),  3 * v(n) / n, v(4 * a * n), etc.,

        if not expr.has(AppliedUndef):
            self.error()

        const, expr = factor_const(expr, n)

        if isinstance(expr, AppliedUndef):
            return self.func(expr, n, f) * const

        tsym = sympify(str(n))
        expr = expr.subs(tsym, n)

        rest = sym.S.One
        undefs = []
        for factor in expr.as_ordered_factors():
            if isinstance(factor, AppliedUndef):
                if factor.args[0] != n:
                    self.error('Weird function %s not of %s' % (factor, n))
                undefs.append(factor)
            else:
                rest *= factor

        if rest.has(AppliedUndef):
            # Have something like 1/v(n)
            self.error()

        exprs = undefs
        if rest.has(n):
            exprs = exprs + [rest]
            rest = sym.S.One

        result = self.term(exprs[0], n, f) * rest

        if len(exprs) == 1:
            return result * const

        self.error('TODO')

    def term(self, expr, n, f):

        const, expr = factor_const(expr, n)
        args = expr.args
        twopidt = 2 * sym.pi * dt
        xn_fac = []
        
        # Check for constant.
        if not expr.has(n):
            return self.add_images(expr * DiracDelta(f) * const, f) / dt

        elif expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            return self.function(expr, n, f) * const
        
        # Handle step u(n-n0)   
        elif expr.is_Function and expr.func in (sym.Heaviside, UnitStep):
            if args[0] is n:
                return const / (1 - sym.exp(-sym.I * twopidt * f)) + const * self.add_images(DiracDelta(f), f) / dt / 2   
            else:
                delay = n - args[0]
                if not delay.has(n):
                    return  const * sym.exp(-sym.I * delay * twopidt * f) / (1 - sym.exp(-sym.I * twopidt * f)) + const * self.add_images(DiracDelta(f), f) / dt / 2
        
        # Handle impulse delta (n-n0)   
        elif expr.is_Function and expr.func == UnitImpulse:
            if args[0] is n:
                return const   
            else:
                delay = n - args[0]
                if not delay.has(n):
                    return  const * sym.exp(-sym.I * delay * twopidt * f)

        # Handle signum
        elif (len(args) == 1 and expr.is_Function and
              expr.func == sym.sign and (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)  
            delay = bb / aa
            if delay.is_integer:
                return const * 2 * sym.exp(sym.I * twopidt * f * delay) / (1 - sym.exp(-sym.I * twopidt * f))
    
    
        # Handle sincu        
        elif (len(args) == 1 and expr.is_Function and
              (expr.func == sincu or expr.func == sincn) and
              (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1) 
            bb = args[0].coeff(n, 0) 
            delay = bb/aa
            fac1 = sym.pi / aa  
            fac2 = 1 /  twopidt
            if expr.func == sincn:
                fac1 = 1 / aa
                fac2 = sym.pi / twopidt
            if delay.is_integer:
                return const * fac1 * (UnitStep(f + fac2 * aa)  - UnitStep(f - fac2 * aa) ) *sym.exp(sym.I * delay * twopidt * f)
        
        # Handle rect
        elif (len(args) == 1 and expr.is_Function and expr.func == rect and
            (args[0].as_poly(n)).is_linear):
            qq = 1/ args[0].coeff(n, 1)
            delay = qq * args[0].coeff(n, 0)
            qq = qq // 2
            if delay.is_integer:
                return const * sym.exp(sym.I * delay * twopidt * f) * sym.sin(twopidt * f / 2 * (2 * qq + 1) ) / sym.sin(twopidt * f / 2)
        
        # Handle cos(a*n+b) 
        elif (len(args) == 1 and expr.is_Function and
            expr.func == sym.cos and (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1) / twopidt
            bb = args[0].coeff(n, 0)
            ret = sym.exp(-sym.I * bb) * DiracDelta(f + aa) + sym.exp(sym.I * bb) * DiracDelta(f - aa)
            return const * self.add_images(ret, f) / dt / 2
        
        # Handle sin(a*n+b) 
        elif (len(args) == 1 and expr.is_Function
            and expr.func == sym.sin and (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1) / twopidt
            bb = args[0].coeff(n, 0) 
            ret = sym.exp(-sym.I * bb) * DiracDelta(f + aa) - sym.exp(sym.I * bb) * DiracDelta(f - aa)
            return const * sym.I * self.add_images(ret, f) / dt / 2

        # Handle exp(j*a*n+b)
        elif (len(args) == 1 and expr.is_Function and 
                expr.func == sym.exp and (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1) / twopidt / sym.I
            bb = args[0].coeff(n, 0) 
            # Check if for complex exp function with abs value to
            # allow for parameter arguments
            if abs(expr / sym.exp(bb)) == 1:
                return const * sym.exp(bb) * self.add_images(DiracDelta(f - aa), f) / dt         
            
        # Multiplication with n       use n* x(n)  o--o  j / twopidt * d/df X(f)
        elif is_multiplied_with(expr, n, 'n', xn_fac):
            expr = expr / xn_fac[0]
            X = self.transform(expr, n, f)
            return const / twopidt * sym.I * sym.simplify(sym.diff(X, f))

        # Multiplication with a**(b*n+c)   use    lam**n * x(n)  o--o  X(f/lam)
        elif is_multiplied_with(expr, n, 'a**n', xn_fac):
            expr /= xn_fac[0]
            ref = xn_fac[0].args
            lam = ref[0]
            bb = ref[1].coeff(n, 1)
            cc = ref[1].coeff(n, 0)              
            X = self.term(expr, n, f)
            return lam**cc * sym.simplify(X.subs(f, f / lam**bb))

        # Multiplication with exp(b*n+c)   use    exp**n * x(n)  o--o  X(f/exp(1))
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac):
            expr /= xn_fac[0]
            ref = xn_fac[0].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)              
            X = self.term(expr, n, f)
            return sym.exp(cc) * sym.simplify(X.subs(f, f / sym.exp(bb))) 
        
        return const * self.sympy(expr, n, f)

    
dtft_transformer = DTFTTransformer()


def discrete_time_fourier_transform(expr, n, f, images=0, evaluate=True,
                                    **assumptions):
    """Compute bilateral discrete-time Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(f)
    """

    return dtft_transformer.transform(expr, n, f, evaluate=evaluate,
                                      images=images, **assumptions)


def DTFT(expr, n, f, images=0, evaluate=True, **assumptions):
    """Compute bilateral discrete-time Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(f)
    """

    return dtft_transformer.transform(expr, n, f, evaluate=evaluate,
                                      images=images, **assumptions)    

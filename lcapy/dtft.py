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
        msym = symsymbol('m', integer=True)        
        foo = expr.replace(f, f - msym / dt)
        result = sym.summation(foo, (msym, self.m1, self.m2))
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
        Omega = 2 * sym.pi * f * dt
        
        # Check for constant.
        if not expr.has(n):
            return self.add_images(expr * DiracDelta(f) * const, f)

        if expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            return self.function(expr, n, f) * const

        # handle cos(a*n+b) 
        if (len(args) == 1 and expr.is_Function and
            expr.func == sym.cos and (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)
            om_0 = expr.args
            co = sym.cos(bb) * (DiracDelta(Omega + aa) + DiracDelta(Omega - aa))
            si = -sym.sin(bb) * (DiracDelta(Omega + aa) - DiracDelta(Omega - aa))
            return self.add_images(sym.pi * (co + sym.I * si) * const, f)
        
        # handle sin(a*n+b) 
        if (len(args) == 1 and expr.is_Function
            and expr.func == sym.sin and (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)
            om_0 = expr.args
            co = sym.sin(bb) * (DiracDelta(Omega + aa) + DiracDelta(Omega - aa))
            si = sym.cos(bb) * (DiracDelta(Omega + aa) - DiracDelta(Omega - aa))
            return self.add_images(sym.pi * (co + sym.I * si) * const, f)
        
        # handle signum
        if (len(args) == 1 and expr.is_Function and
            expr.func == sym.sign and (args[0].as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)  
            delay = bb / aa
            if delay.is_integer:
                return 2 * sym.exp(sym.I * Omega * delay) / (1 - sym.exp(-sym.I * Omega))

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

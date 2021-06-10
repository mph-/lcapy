"""This module provides support for discrete-time Fourier transforms. 

 It calculates the discrete-times Fourier transform using:

   X(f) = \sum_{n=-M}^{M} x(n) e^{-j * 2 * \pi * n * dt * f}

where 2 * M is the number of images.

Copyright 2021 Michael Hayes, UCECE

"""

import sympy as sym
from sympy import oo, DiracDelta
from .transformer import Transformer
from .sym import sympify, AppliedUndef, j, pi
from .dsym import dt
from .extrafunctions import UnitImpulse, UnitStep
from .utils import factor_const, scale_shift
from .matrix import Matrix

__all__ = ('DTFT', )


class DTFTTransformer(Transformer):

    name = 'DTFT'
    inverse = False
    
    def key(self, expr, n, f, **assumptions):
        return expr, n, f

    def noevaluate(self, expr, n, f):

        foo = expr * sym.exp(-2 * j * pi * n * dt * f)
        result = sym.Sum(foo, (n, -oo, oo))
        return result

    def check(self, expr, n, f, images, **assumptions):

        # Hack
        self.images = images
        self.m1 = -(images // 2)
        self.m2 = self.m1 + self.images
        
        if expr.has(f):
            self.error('Expression depends on f')
        
        if expr.is_Piecewise and expr.args[0].args[1].has(n >= 0):
            self.error('Expression is unknown for n < 0' % expr)
    
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
        if self.inverse:
            func = name[0].lower() + name[1:] + '(%s)' % f
        else:
            func = name[0].upper() + name[1:] + '(%s)' % f

        result = sympify(func).subs(fsym, f / scale) / abs(scale)

        if shift != 0:
            if self.inverse:
                shift = -shift
            result = result * sym.exp(2 * sym.I * sym.pi * f * shift / scale)

        return result

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

    def term1(self, expr, n, f):

        const, expr = factor_const(expr, n)

        # Check for constant.
        if not expr.has(n):
            # TODO, add images.
            return expr * DiracDelta(f) * const

        if expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            return self.function(expr, n, f) * const

        return const * self.sympy(expr, n, f)

    def term(self, expr, n, f):

        result = self.term1(expr, n, f)
        if self.inverse:
            result /= self.N
        return result
    
    
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

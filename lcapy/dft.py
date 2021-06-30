"""This module provides support for discrete Fourier transforms. 

 It calculates the discrete Fourier transform using:

   X(k) = \sum_{n=0}^{N-1} x(n) e^{-j * 2 * \pi * n * k / N}

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy as sym
from .transformer import BilateralForwardTransformer
from .sym import sympify, AppliedUndef, j, pi
from .extrafunctions import UnitImpulse, UnitStep
from .utils import factor_const, scale_shift
from .matrix import Matrix

__all__ = ('DFT', 'DFTmatrix')


class DFTTransformer(BilateralForwardTransformer):

    name = 'DFT'
    is_inverse = False
    
    def key(self, expr, n, k, **assumptions):
        return expr, n, k, assumptions.get('N', None)

    def noevaluate(self, expr, n, k):

        foo = expr * sym.exp(-2 * j * pi * n * k / self.N)
        result = sym.Sum(foo, (n, 0, self.N - 1))
        return result

    def check(self, expr, n, k, N=None, **assumptions):

        self.N = N
        
        if expr.has(k):
            self.error('Expression depends on k')
        
        if expr.is_Piecewise and expr.args[0].args[1].has(n >= 0):
            self.error('Expression is unknown for n < 0 (use causal=True)')
    
    def sympy(self, expr, n, k):

        foo = expr * sym.exp(-2 * j * pi * n * k / self.N)
        result = sym.summation(foo, (n, 0, self.N - 1))

        return result

    def func(self, expr, n, k):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], n)

        fsym = sympify(str(k))

        # Convert v(n) to V(k), etc.
        name = expr.func.__name__
        if self.is_inverse:
            func = name[0].lower() + name[1:] + '(%s)' % k
        else:
            func = name[0].upper() + name[1:] + '(%s)' % k

        result = sympify(func).subs(fsym, k / scale) / abs(scale)

        if shift != 0:
            if self.is_inverse:
                shift = -shift
            result = result * sym.exp(2 * sym.I * sym.pi * k * shift / scale)

        if self.is_inverse:
            result *= self.N

        return result


    def function(self, expr, n, k):

        # Handle expressions with a function of FOO, e.g.,
        # v(n), v(n) * y(n),  3 * v(n) / n, v(4 * a * n), etc.,

        if not expr.has(AppliedUndef):
            self.error()

        const, expr = factor_const(expr, n)

        if isinstance(expr, AppliedUndef):
            return self.func(expr, n, k) * const

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

        result = self.term(exprs[0], n, k) * rest

        if len(exprs) == 1:
            return result * const

        dummy = 'm' if self.is_inverse else 'l'

        for m in range(len(exprs) - 1):
            if m == 0:
                nu = sympify(dummy)
            else:
                nu = sympify(dummy + '_%d' % m)
            expr2 = self.term(exprs[m + 1], n, k)
            # Should be a circular convolution.
            result = sym.Sum(result.subs(k, k - nu) * expr2.subs(k, nu),
                             (nu, 0, self.N - 1)) / self.N

        return result * const

    def term1(self, expr, n, k):

        const, expr = factor_const(expr, n)

        # Check for constant.
        if not expr.has(n):
            return expr * self.N * UnitImpulse(k) * const

        if expr.is_Function and expr.func == UnitStep and expr.args[0] == n:
            return const * self.N * UnitImpulse(k)    

        if expr.is_Function and expr.func == UnitImpulse and expr.has(n):
            scale, shift = scale_shift(expr.args[0], n)
            # This implies a circular shift rather than a delay/advance.
            return const * sym.exp(2 * j * pi * shift * k / self.N)

        if expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            return self.function(expr, n, k) * const

        if expr.is_Function and expr.func == sym.exp and expr.args[0].has(n):
            scale, shift = scale_shift(expr.args[0], n)
            l = scale * self.N / (j * 2 * pi)
            # TODO, need to wrap the frequency shift, l
            return const * self.N * sym.exp(2 * j * pi * shift * k / self.N) * UnitImpulse(k + l)

        return const * self.sympy(expr, n, k)

    def term(self, expr, n, k):

        result = self.term1(expr, n, k)
        if self.is_inverse:
            result /= self.N
        return result
    
    
dft_transformer = DFTTransformer()


def discrete_fourier_transform(expr, n, k, N=None, evaluate=True,
                               **assumptions):
    """Compute bilateral discrete Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(k)
    """

    return dft_transformer.transform(expr, n, k, evaluate=evaluate, N=N,
                                     **assumptions)


def DFT(expr, n, k, N=None, evaluate=True, **assumptions):
    """Compute bilateral discrete Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(k)
    """

    return dft_transformer.transform(expr, n, k, evaluate=evaluate, N=N,
                                     **assumptions)    


def DFTmatrix(N):
    """Return DFT matrix of size `N` x `N`."""    

    from .functions import exp
    from .sym import j, pi
    
    w = exp(-j * 2 * pi / N)

    a = Matrix.zeros(N)
    
    for row in range(N):
        for col in range(N):
            a[row, col] = w ** (row * col)
    return a

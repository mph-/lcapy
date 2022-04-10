"""This module provides support for the inverse DFT.  It
calculates the bilateral inverse DFT transform using:

   s(n) = (1 / N) * \sum_{n=0}^{N-1} S(k) e^{j * 2 * \pi * n * k / N}

Copyright 2021 Michael Hayes, UCECE

"""

from .dft import DFTTransformer
from .matrix import Matrix
from .sym import j, pi
import sympy as sym


__all__ = ('IDFT', 'IDFTmatrix')


class InverseDFTTransformer(DFTTransformer):

    name = 'inverse DFT'
    is_inverse = True

    def noevaluate(self, expr, k, n):

        foo = expr * sym.exp(2 * j * pi * n * k / self.N)
        result = sym.Sum(foo, (n, 0, self.N - 1))
        result /= self.N
        return result

    def check(self, expr, k, n, N=None, **assumptions):

        if expr.has(n):
            self.error('Expression depends on n')
        super(InverseDFTTransformer, self).check(expr, k, n, N, **assumptions)


inverse_dft_transformer = InverseDFTTransformer()


def inverse_discrete_fourier_transform(expr, k, n, N, evaluate=True, **assumptions):
    """Compute bilateral inverse DFT transform of expr.

    Undefined functions such as X(k) are converted to x(n)
    """

    return inverse_dft_transformer.transform(expr, k, n, evaluate=evaluate,
                                             N=N, **assumptions)


def IDFT(expr, k, n, N, evaluate=True, **assumptions):
    """Compute bilateral inverse discrete Fourier transform of expr.

    Undefined functions such as X(k) are converted to x(n)
    """

    return inverse_dft_transformer.transform(expr, k, n, evaluate=evaluate,
                                             N=N, **assumptions)


def IDFTmatrix(N):
    """Return inverse DFT matrix of size `N` x `N`."""

    from .functions import exp

    w = exp(j * 2 * pi / N)

    a = Matrix.zeros(N)

    for row in range(N):
        for col in range(N):
            a[row, col] = w ** (row * col)
    return a / N

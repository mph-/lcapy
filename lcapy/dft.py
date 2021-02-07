"""This module provides support for discrete Fourier transforms. 

 It calculates the discrete Fourier transform using:

   X(k) = \sum_{n=0}^{N-1} x(n) e^{-j * 2 * \pi * n * k / N}

Copyright 2020 Michael Hayes, UCECE

"""

import sympy as sym
from .sym import sympify, AppliedUndef, j, pi
from .functions import UnitImpulse, UnitStep, exp
from .utils import factor_const, scale_shift
from .matrix import Matrix

__all__ = ('DFT', 'IDFT', 'DFTmatrix', 'IDFTmatrix')


discrete_fourier_cache = {}


def discrete_fourier_sympy(expr, n, k, N):

    foo = expr * sym.exp(-2 * j * pi * n * k / N)
    result = sym.summation(foo, (n, 0, N - 1))
    
    return result


def discrete_fourier_func(expr, n, k, N, inverse=False):

    if not isinstance(expr, AppliedUndef):
        raise ValueError('Expecting function for %s' % expr)

    scale, shift = scale_shift(expr.args[0], n)

    fsym = sympify(str(k))
    
    # Convert v(n) to V(k), etc.
    name = expr.func.__name__
    if inverse:
        func = name[0].lower() + name[1:] + '(%s)' % k
    else:
        func = name[0].upper() + name[1:] + '(%s)' % k

    result = sympify(func).subs(fsym, k / scale) / abs(scale)

    if shift != 0:
        if inverse:
            shift = -shift
        result = result * sym.exp(2 * sym.I * sym.pi * k * shift / scale)

    if inverse:
        result *= N
        
    return result


def discrete_fourier_function(expr, n, k, N, inverse=False):

    # Handle expressions with a function of FOO, e.g.,
    # v(n), v(n) * y(n),  3 * v(n) / n, v(4 * a * n), etc.,
    
    if not expr.has(AppliedUndef):
        raise ValueError('Could not compute Fourier transform for ' + str(expr))

    const, expr = factor_const(expr, n)
    
    if isinstance(expr, AppliedUndef):
        return discrete_fourier_func(expr, n, k, N, inverse) * const
    
    tsym = sympify(str(n))
    expr = expr.subs(tsym, n)

    rest = sym.S.One
    undefs = []
    for factor in expr.as_ordered_factors():
        if isinstance(factor, AppliedUndef):
            if factor.args[0] != n:
                raise ValueError('Weird function %s not of %s' % (factor, n))
            undefs.append(factor)
        else:
            rest *= factor

    if rest.has(AppliedUndef):
        # Have something like 1/v(n)
        raise ValueError('Cannot compute Fourier transform of %s' % rest)
            
    exprs = undefs
    if rest.has(n):
        exprs = exprs + [rest]
        rest = sym.S.One

    result = discrete_fourier_term(exprs[0], n, k, N, inverse) * rest
        
    if len(exprs) == 1:
        return result * const

    dummy = 'm' if inverse else 'l'

    for m in range(len(exprs) - 1):
        if m == 0:
            nu = sympify(dummy)
        else:
            nu = sympify(dummy + '_%d' % m)
        expr2 = discrete_fourier_term(exprs[m + 1], n, k, N, inverse)
        # Should be a circular convolution.
        result = sym.Sum(result.subs(k, k - nu) * expr2.subs(k, nu),
                         (nu, 0, N - 1)) / N
    
    return result * const

def discrete_fourier_term(expr, n, k, N, inverse=False):

    const, expr = factor_const(expr, n)

    # Check for constant.
    if not expr.has(n):
        return expr * N * UnitImpulse(k) * const

    if expr.is_Function and expr.func == UnitStep and expr.args[0] == n:
        return const * N * UnitImpulse(k)    

    if expr.is_Function and expr.func == UnitImpulse and expr.has(n):
        scale, shift = scale_shift(expr.args[0], n)
        # This implies a circular shift rather than a delay/advance.
        return const * sym.exp(2 * j * pi * shift * k / N)

    if expr.has(AppliedUndef):
        # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
        return discrete_fourier_function(expr, n, k, N, inverse) * const

    if expr.is_Function and expr.func == sym.exp and expr.args[0].has(n):
        scale, shift = scale_shift(expr.args[0], n)
        l = scale * N / (j * 2 * pi)
        # TODO, need to wrap the frequency shift, l
        return const * N * sym.exp(2 * j * pi * shift * k / N) * UnitImpulse(k + l)
    
    return const * discrete_fourier_sympy(expr, n, k, N)


def discrete_fourier_transform(expr, n, k, N=None, inverse=False, evaluate=True):
    """Compute bilateral discrete Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(k)
    """

    if expr.is_Equality:
        return sym.Eq(discrete_fourier_transform(expr.args[0], n, k, N, inverse),
                      discrete_fourier_transform(expr.args[1], n, k, N, inverse))    

    if not evaluate:
        foo = expr * sym.exp(-2 * j * pi * n * k / N)
        result = sym.Sum(foo, (n, 0, N - 1))
        if inverse:
            result /= N
        return result
    
    key = (expr, n, k, N, inverse)
    if key in discrete_fourier_cache:
        return discrete_fourier_cache[key]

    if not inverse and expr.has(k):
        raise ValueError('Cannot discrete Fourier transform for expression %s that depends on %s' % (expr, k))


    if expr.is_Piecewise and expr.args[0].args[1].has(n >= 0):
        raise ValueError('Cannot discrete Fourier transform for expression %s that is unknown for n < 0' % expr)
    
    if inverse:
        n, k = k, n

    var = sym.Symbol(str(n))
    expr = expr.replace(var, n)

    orig_expr = expr

    if expr.has(sym.cos) or expr.has(sym.sin):
        expr = expr.rewrite(sym.exp)
        
    terms = expr.expand().as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += discrete_fourier_term(term, n, k, N, inverse=inverse)
    except ValueError:
        raise ValueError('Could not compute discrete Fourier transform for ' + str(orig_expr))

    if inverse:
        result /= N
    
    discrete_fourier_cache[key] = result
    return result


def inverse_discrete_fourier_transform(expr, k, n, N=None, evaluate=False):
    """Compute bilateral inverse discrete Fourier transform of expr.

    Undefined functions such as X(k) are converted to x(n)
    """
    
    if expr.has(n):
        raise ValueError('Cannot inverse discrete Fourier transform for expression %s that depends on %s' % (expr, n))
    
    result = discrete_fourier_transform(expr, n, k, N, inverse=True, evaluate=evaluate) 
    return sym.simplify(result)


def DFT(expr, n, k, N, **assumptions):
    """Compute bilateral discrete Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(k)
    """
    
    return discrete_fourier_transform(expr, n, k, N, **assumptions)


def IDFT(expr, k, n, N, **assumptions):
    """Compute bilateral inverse discrete Fourier transform of expr.

    Undefined functions such as X(k) are converted to x(n)
    """    
    
    return inverse_discrete_fourier_transform(expr, k, n, N, **assumptions)


def DFTmatrix(N):
    """Return DFT matrix of size `N` x `N`."""    

    w = exp(-j * 2 * pi / N)

    a = Matrix.zeros(N)
    
    for row in range(N):
        for col in range(N):
            a[row, col] = w ** (row * col)
    return a


def IDFTmatrix(N):
    """Return inverse DFT matrix of size `N` x `N`."""

    w = exp(j * 2 * pi / N)

    a = Matrix.zeros(N)
    
    for row in range(N):
        for col in range(N):
            a[row, col] = w ** (row * col)
    return a / N

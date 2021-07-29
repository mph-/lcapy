"""This module provides support for discrete Fourier transforms. 

 It calculates the discrete Fourier transform using:

   X(k) = \sum_{n=0}^{N-1} x(n) e^{-j * 2 * \pi * n * k / N}

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy as sym
from .transformer import BilateralForwardTransformer
from .sym import sympify, AppliedUndef, j, pi, symsymbol
from .extrafunctions import UnitImpulse, UnitStep
from .utils import factor_const, scale_shift
from .matrix import Matrix
from .ztransform import is_multiplied_with

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

        try:
            N = N.expr
        except:
            pass
        
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
                nu = symsymbol(dummy, integer=True)
            else:
                nu = symsymbol(dummy + '_%d' % m, integer=True)
            expr2 = self.term(exprs[m + 1], n, k)
            # Should be a circular convolution.
            result = sym.Sum(result.subs(k, k - nu) * expr2.subs(k, nu),
                             (nu, 0, self.N - 1)) / self.N

        return result * const

    
    # Make transform  Xq = sum_lower^upper  x[n] * q**n
    # Returns Xq and cases.
    # Cases contains special Xq values at e.g. k=0 or k=1 depending on x[n]
    # May also be useful for the Z-transform or the DTFT for finite sequences. 
    def term1(self, expr, n, q, lower, upper):

        const, expr = factor_const(expr, n)    
        args = expr.args
        xn_fac = []    

        # Check for constant.
        if not expr.has(n):
            result_q = const * expr * q**lower * (1 - q**(upper - lower + 1)) / (1 - q)
            # Special case k = 0
            cases = {0 : const * expr * (upper - lower + 1)}
            return result_q, cases       

        # Handle delta(n-n0)
        elif expr.is_Function and expr.func == UnitImpulse and ((expr.args[0]).as_poly(n)).is_linear:
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0) 
            nn0 = -bb / aa
            # FIXME since upper depends on unknown N
            # if nn0 > upper or nn0 < lower:
            #     return 0 * q, {}
            if nn0.is_integer or nn0.is_integer is None:
                return_q = const * q**nn0
                return return_q, {}             

        # Handle n**p
        if (expr == n or
            (expr.is_Pow and args[1].is_integer and args[1].is_positive
             and args[0] == n)):
            p = 1
            try:
                p = args[1]
            except:
                pass
            # derivatives          
            result_q = (q**lower - q**(upper + 1)) / (1 - q) 
            for i in range(p):
                result_q = q * sym.diff(result_q, q)             
            result_q *= const

            # Special case k=0, use Faulhaber's formula
            cases = {0 : const * sym.factor((sym.bernoulli(p + 1, upper + 1) - sym.bernoulli(p + 1, lower)) / (p + 1))}            
            return result_q, cases      

        # Handle  *rect((n-a)/b)
        elif is_multiplied_with(expr, n, 'rect', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = 1 / sym.expand(ref[0]).coeff(n, 1) 
            aa = -bb * sym.expand(ref[0]).coeff(n, 0)
            # left and right  index
            nn0 = aa - bb // 2
            nn1 = nn0 + bb - 1            
            if nn0.is_integer:
                if nn0 > upper:
                    return 0 * q, {}
                else:
                    nn0 = max(lower, nn0)
            if nn1.is_integer:
                if nn1 < lower:                
                    return 0 * q, {}
                else:
                    nn1 = min(upper, nn1)
            Xq, cases = self.term1(expr, n, q, nn0, nn1)
            return_q = const * Xq
            for key in cases:
                cases[key] = const * cases[key]             
            return return_q, cases

        # Handle  *u(n-n0)
        elif is_multiplied_with(expr, n, 'UnitStep', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            aa = ref[0].coeff(n, 1) 
            bb = ref[0].coeff(n, 0)
            if abs(aa) != 1:
                print("Use u(n-n0)")
            else:    
                # Positive step
                if aa == 1 and ((bb.is_integer and -bb < self.N and -bb >= 0) or not bb.is_number):
                    Xq, cases = self.term1(expr, n, q, -bb, self.N - 1)
                # Negative step    
                elif aa == -1 and ((bb.is_integer and bb < self.N and bb > 0) or not bb.is_number):
                    Xq, cases = self.term1(expr, n, q, 0, bb)
                else:
                    Xq = 0 * q
                    cases = {}

                result_q = const * Xq
                for key in cases:
                    cases[key] = const * cases[key]                  

                return result_q, cases    

        # Handle  *exp(j*a*n+b) 
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac) and abs(xn_fac[-1] / sym.exp(args[0].coeff(n, 0))) == 1:
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            aa = sym.expand(ref[0]).coeff(n, 1) / sym.I
            bb = sym.expand(ref[0]).coeff(n, 0) 
            # find transform
            Xq, ca =  self.term1(expr, n, q, lower, upper)
            # check frequency
            if aa.is_constant() and abs(aa) >= pi:
                print("Warning: Frequency may be out of range")                
            result_q = const * sym.exp(bb) * Xq.subs(q, q * sym.exp(sym.I * aa))
            cases = {}
            # check special case and shift accordingly
            k0 = aa * self.N / 2 / pi
            if k0.is_integer and len(ca) != 0:
                for key in ca:
                    cases[key + k0] = const * ca[key] * sym.exp(bb)
            return result_q, cases         

        # Handle *sin(b*n+c)
        elif is_multiplied_with(expr, n, 'sin(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0) 
            # check frequency
            if abs(bb) >= pi:
                print("WARNING: Frequency may be out of range")             
            Xq, ca =  self.term1(expr, n, q, lower, upper)
            Xq1 = sym.exp(sym.I * cc) * Xq.subs(q, q * sym.exp(sym.I * bb))
            Xq2 = sym.exp(-sym.I * cc) * Xq.subs(q, q * sym.exp(-sym.I * bb))
            result_q = const * (Xq1 - Xq2) / 2 /sym.I
            cases = {}
            # Check special case shift abs
            k0 = bb * self.N / 2 / pi
            if k0.is_integer and len(ca) != 0:
                # Handle special cases and shift left and right
                for key in ca:
                    cases[key + k0] = const * ca[key] * sym.exp(sym.I * cc) / 2 / sym.I                     
                    cases[key - k0] = -const * ca[key] * sym.exp(-sym.I * cc) / 2 / sym.I
            return  result_q, cases

        # Handle *cos(b*n+c)
        elif is_multiplied_with(expr, n, 'cos(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0) 
            # Check frequency
            if abs(bb) >= pi:
                print("Warning: Frequency may be out of range")             
            Xq, ca =  self.term1(expr, n, q, lower, upper)
            Xq1 = sym.exp(sym.I * cc) * Xq.subs(q, q * sym.exp(sym.I * bb))
            Xq2 = sym.exp(-sym.I * cc) * Xq.subs(q, q * sym.exp(-sym.I * bb))
            result_q = const * (Xq1 + Xq2) / 2
            cases = {}
            # Handle special cases if necessary and shift left and right
            k0 = bb * self.N / 2 / pi
            if k0.is_integer and len(ca) != 0:
                for key in ca:
                    cases[key + k0] = const * ca[key] * sym.exp(sym.I * cc) / 2 
                    cases[key - k0] = const * ca[key] * sym.exp(-sym.I * cc) / 2
            return result_q, cases        

        # Handle  *exp(bb*n+cc)
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac):
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)                 
            Xq, ca =  self.term1(expr, n, q, lower, upper)
            result_q = const * sym.exp(cc) * Xq.subs(q, q * sym.exp(bb))                     
            # No special cases 
            cases = {}
            return result_q, cases

        # Handle  *a**n
        elif is_multiplied_with(expr, n, 'a**n', xn_fac):
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            lam = ref[0]
            bb = ref[1].coeff(n, 1)
            cc = ref[1].coeff(n, 0) 
            Xq, Xk =  self.term1(expr, n, q, lower, upper)
            result_q = const * lam**cc * Xq.subs(q, q * lam ** bb)
            # No special cases 
            cases = {}
            return result_q, cases        

        # Handle *n       
        elif is_multiplied_with(expr, n, 'n', xn_fac):
            expr = expr / xn_fac[-1]
            Xq, ca = self.term1(expr, n, q, lower, upper)
            result_q = const * q * sym.diff(Xq, q) 
            cases = {}
            if len(ca) != 0:
                raise ValueError("No sym.diff possible for discrete cases")
            return result_q, cases                

        return const * sym.summation(expr * q**n, (n, lower, upper)), {}    
    
    def term(self, expr, n, k):

        if expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            result = self.function(expr, n, k)
            if self.is_inverse:
                result /= self.N
            return result
        
        q = sym.Symbol('q')

        Xq, ca = self.term1(expr, n, q, 0, self.N - 1)
        result = Xq.subs(q, sym.exp(-sym.I * 2 * pi / self.N * k))
        # Add special cases with delta(k-k0)*values
        for key in ca:
            k0 = key
            if k0 < 0:
                # Shift negative frequencies
                k0 += self.N 
            result += UnitImpulse(k-k0) * ca[key].subs(q, sym.exp(-sym.I * 2 * pi / self.N * k0))
        
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

"""This module provides support for discrete Fourier transforms. 

It calculates the discrete Fourier transform using:

   X(k) = \sum_{n=0}^{N-1} x(n) e^{-j * 2 * \pi * n * k / N}

Most of the cunningness in this file is from Juergen Weizenecker.

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy as sym
from .transformer import BilateralForwardTransformer
from .sym import sympify, AppliedUndef, j, pi, symsymbol
from .extrafunctions import UnitImpulse, UnitStep
from .utils import factor_const, scale_shift
from .matrix import Matrix
from .ztransform import is_multiplied_with
from copy import deepcopy

__all__ = ('DFT', 'DFTmatrix')


def is_in_interval(n0, n1, nk, comment=''):
    """Check if bin belongs to an interval [n0, n1]. 
     n0, n1, nk may contain parameter
        
    This returns:
     -1  if bin is not in interval and is on the left side
      0  if bin is in interval
      1  if bin is not in interval and is on the right side

    If n0, n1 or nk contain parameter, 0 is returned if no decision is possible and
    a comment is printed.
    """
    
    # Is in [n0, n1]
    if sym.sympify(nk - n0).is_nonnegative and sym.sympify(n1 - nk).is_nonnegative:
        return 0
    # Is on the left
    elif sym.sympify(nk - n0).is_negative:
        return -1
    # Is on the right
    elif sym.sympify(n1 - nk).is_negative:
        return 1
    # No decision, assume in interval
    else:
        print(comment, 'assume bin', nk, 'is in interval [%s, %s]' % (n0, n1))
        return 0 

    
def find_pow(q_expr, base, NN):
    """This finds terms of the form base**NN in the expression q_expr.
    It is used for simplifying q**N terms to one if q=exp(-+j*2*pi/N*k).
    
    It returns a list of all bases found e.g., [q/3, q/5, 2q, aq,  .......].
    """
    
    result = []
    for x in q_expr.atoms(sym.Pow):
        this_ex = x.as_base_exp()[1]
        this_ba = x.as_base_exp()[0]
        if this_ex == NN and not sym.simplify(this_ba / base).has(base):
            result += [this_ba] 
    
    return result


class QkTransform(object):
    """This handles X(q) transforms, where

    X(q) = sum_lower^upper x(n) * q**N  

    As special cases are possible (e.q. for q=0) the class handles the
    general expression and special cases"""
    
    def __init__(self, val_q, k0=None, val0=None, valqk=None, val_k=0):
        self.Xq = val_q
        self.Xk = val_k
        if k0 is not None:
            self.has_special = True
            self.cases = [k0]
            self.case_expr = {k0: [val0, valqk]} 
        else:
            self.has_special = False
            self.cases = []
            self.case_expr = {}             
            
    def multiply(self, c0):
        """Multiply by constant."""
        
        self.Xq *= c0
        self.Xk *= c0
        if self.has_special:
            for key in self.case_expr:
                self.case_expr[key][0] *= c0
                self.case_expr[key][1] *= c0
     
    def rm_cases(self):
        """Delete special cases."""

        if self.has_special:
            self.has_special = False
            self.cases = []
            self.case_expr = {}            
    
    def shift_k(self, dk):
        """Shift the special cases k values by dk."""
        
        if self.has_special:
            self.cases = [ii+dk for ii in self.cases]
            new_dict = {}
            for key in self.case_expr:
                new_dict[key+dk] = self.case_expr[key]
            self.case_expr = new_dict 
    
    def add(self, sec):
        """Merge (add) two QkTransforms."""
        
        self.Xq += sec.Xq
        self.Xk += sec.Xk
        # Both have special cases
        if sec.has_special and self.has_special:            
            self.cases += [ii for ii in sec.cases if ii not in self.cases]
            for key in sec.case_expr:
                if key in self.case_expr:
                    self.case_expr[key][0] += sec.case_expr[key][0]
                    self.case_expr[key][1] += sec.case_expr[key][1]
                else:
                    self.case_expr[key]= sec.case_expr[key]
        # Only the second has special cases  
        elif sec.has_special and  not self.has_special:
            self.has_special = True
            self.cases = [ii for ii in sec.cases]
            for key in sec.case_expr:
                self.case_expr[key] = sec.case_expr[key]
    
    def subs(self, q, q_new):
        """Substitute q by a new expression."""
        
        self.Xq = (self.Xq).subs(q, q_new)
        if self.has_special:
            for key in self.case_expr:
                self.case_expr[key][1] = self.case_expr[key][1].subs(q, q_new)        
    

    def simp_qN(self, ba, NN, all=False):
        """Simplifies ba**NN to 1 
        
        If all is False, only the special cases will regarded.
        Replace e.g., (q/3)**NN by (1/3)**N    
        """
        
        # Simplify for "no special cases" or if demanded
        if all or not self.has_special:
            # General expresion        
            terms = find_pow(self.Xq, ba, NN)            
            for full_base in terms:
                self.Xq = self.Xq.replace(full_base**NN, (full_base / ba)**NN)
            
        # Special cases    
        if self.has_special:
            for key in self.case_expr:
                terms = find_pow(self.case_expr[key][1], ba, NN)
                for full_base in terms:
                    self.case_expr[key][1] = self.case_expr[key][1].replace(full_base**NN, (full_base / ba)**NN)                
        
            
    def make_transform(self, NN, backward, q, sub, k):
        """Assemble the final tranform from all the cases, the general case, or Xk.
         
        backward is True for IDFT
        q is replaced by sub 
        """

        if self.has_special:
            for ca in self.cases:
                k0 = ca
                # Shift k0 to [0, N-1]
                if backward and sym.sympify(ca).is_positive:
                    k0 -= NN
                elif not backward and sym.sympify(ca).is_negative:
                    k0 += NN
                # Special case * delta 
                self.Xk += self.case_expr[ca][0] * UnitImpulse(k - k0) 
                # Rest * (1-delta) to avoid unwanted superposition
                self.Xk += sym.simplify(self.case_expr[ca][1].subs(q, sub)) * (1 - UnitImpulse(k - k0)) 
        else:
            self.Xk = self.Xq.subs(q, sub) + self.Xk        


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

        # Convert constants to SymPy expressions.
        N = sym.sympify(N)

        if not N.is_integer:
            raise ValueError('%s not integer, redefine with integer=True' % N)
        if not N.is_positive:
            raise ValueError('%s not positive, redefine with positive=True' % N)
        
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
            result = result * sym.exp(2 * sym.I * sym.pi * k * shift / (scale * self.N))

        if self.is_inverse:
            result *= self.N

        return result

    def function(self, expr, n, k):

        # Handle expressions with a function of FOO, e.g., 
        # v(n), v(n) * y(n), 3 * v(n) / n, v(4 * a * n), etc., 

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
    def termXq(self, expr, n, k, q, lower, upper):

        const, expr = factor_const(expr, n)    
        args = expr.args
        xn_fac = []    

        # Check for constant.
        if not expr.has(n):
            result_q = const * expr * q**lower * (1 - q**(upper - lower + 1)) / (1 - q)
            # Special case k = 0
            result_1 = const * expr * (upper - lower + 1)
            result = QkTransform(result_q, 0, result_1, result_q)
            return result       

        # Handle delta(n-n0)
        elif (expr.is_Function and expr.func == UnitImpulse and
              ((expr.args[0]).as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0) 
            nn0 = -bb / aa
            if nn0.is_integer:
                if is_in_interval(lower, upper, nn0, comment='Delta:') == 0:
                    # Shift frequency to -pi/2 ...  pi/2 
                    if nn0.has(self.N):
                        nn0 = nn0.subs(self.N, 0)
                    elif sym.sympify(nn0 - self.N / 2).is_positive:
                        nn0 -= self.N
                        
                    result_q = const * q**nn0
                    result = QkTransform(result_q)
                else:
                    result = QkTransform(0 * q)
                return result
            else:
                print("delta: bin ", nn0, " is not an integer")

        # Handle n**p
        if (expr == n or
            (expr.is_Pow and args[1].is_integer and args[1].is_positive
             and args[0] == n)):
            p = 1
            try:
                p = args[1]
            except:
                pass
            # Derivatives          
            result_q = (q**lower - q**(upper + 1)) / (1 - q) 
            for i in range(p):
                result_q = q * sym.diff(result_q, q)             
            result_q *= const

            # Special case k=0, use Faulhaber's formula
            result_1 = const * sym.factor((sym.bernoulli(p + 1, upper + 1) - sym.bernoulli(p + 1, lower)) / (p + 1))
            result = QkTransform(result_q, 0, result_1, result_q)      
            return result      

        # Handle * rect((n-a)/b)
        elif is_multiplied_with(expr, n, 'rect', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = 1 / sym.expand(ref[0]).coeff(n, 1) 
            aa = -bb * sym.expand(ref[0]).coeff(n, 0)
            # Left and right index
            nn0 = aa - bb // 2
            nn1 = nn0 + bb - 1    
            if (not aa.is_integer) or (not bb.is_integer):
                print("rect((n-a)/b): parameter not an integer")
            else:
                r1 = is_in_interval(lower, upper, nn0, comment='rect left step')
                r2 = is_in_interval(lower, upper, nn1, comment='rect right step')
                if r1 == 0 and r2 == 0:
                    result= self.termXq(expr, n, k, q, nn0, nn1)
                elif r1 == -1 and r2 == 0:
                    result= self.termXq(expr, n, k, q, lower, nn1)
                elif r1 == 0 and r2 == 1:
                    result= self.termXq(expr, n, k, q, nn0, upper)
                elif r1 == -1 and r2 == 1:
                    result= self.termXq(expr, n, k, q, lower, upper)
                else:
                    result = QkTransform(0 * q)
                
                result.multiply(const)             
                return result                 

        # Handle * u(n-n0)
        elif is_multiplied_with(expr, n, 'UnitStep', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            aa = ref[0].coeff(n, 1) 
            bb = ref[0].coeff(n, 0)
            if abs(aa) != 1:
                print("Use u(n-n0)")
            if not bb.is_integer:
                print("Step(n-n0): n0 not an integer")                
            else:    
                # Positive step
                if aa == 1:
                    nn0 = -bb
                    rg = is_in_interval(lower, upper, nn0, comment='Step:')
                    if rg == -1:
                        result= self.termXq(expr, n, k, q, lower, upper)
                    elif rg == 0:
                        result= self.termXq(expr, n, k, q, nn0, upper)
                    else:
                        result = QkTransform(0 * q) 
                        
                # Negative step
                elif aa == -1:
                    nn0 = bb
                    rg = is_in_interval(lower, upper, nn0, comment='Step:')
                    if rg == 1:
                        result= self.termXq(expr, n, k, q, lower, upper)
                    elif rg == 0:
                        result= self.termXq(expr, n, k, q, lower, nn0)
                    else:
                        result = QkTransform(0 * q)

                result.multiply(const)               
                return result  

        # Handle * exp(j*a*n+b) 
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac) and abs(xn_fac[-1] / sym.exp(args[0].coeff(n, 0))) == 1:
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            aa = sym.expand(ref[0]).coeff(n, 1) / sym.I
            bb = sym.expand(ref[0]).coeff(n, 0) 
            # Find transform
            result = self.termXq(expr, n, k, q, lower, upper)
            # Check frequency
            if aa.is_constant() and abs(aa) > pi:
                print("Warning: Frequency may be out of range")       
                
            result.subs(q, q * sym.exp(sym.I * aa))            
            # Check special case and shift accordingly
            k0 = aa * self.N / 2 / pi
            if k0.is_integer and result.has_special:
                result.shift_k(k0)
            else:
                result.rm_cases()
            result.multiply(const * sym.exp(bb))        
            return result         

        # Handle * sin(b*n+c)
        elif is_multiplied_with(expr, n, 'sin(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0) 
            # Check frequency
            if bb.is_constant() and abs(bb) > pi:
                print("Warning: Frequency may be out of range")    
                
            result = self.termXq(expr, n, k, q, lower, upper)
            # Make copies for the transformation:
            # X(q) = (exp(j*bb)*Xq (q*exp(j*bb)) - exp(-j*bb)*Xq(q*exp(-j*bb))) / 2 / j
            rq1 = deepcopy(result)
            rq2 = deepcopy(result)
            # Make general shift
            rq1.subs(q, q * sym.exp(sym.I * bb))           
            rq2.subs(q, q * sym.exp(-sym.I * bb))
            # Check special case shift s
            k0 = bb * self.N / 2 / pi
            if k0.is_integer and result.has_special:
                rq1.shift_k(k0)
                rq2.shift_k(-k0)
            else:
                rq1.rm_cases()
                rq2.rm_cases()
            
            rq1.multiply(const * sym.exp(sym.I * cc) / 2 / sym.I)
            rq2.multiply(-const * sym.exp(-sym.I * cc) / 2 / sym.I)
            
            # Add both parts 
            rq1.add(rq2)
            return  rq1

        # Handle * cos(b*n+c)
        elif is_multiplied_with(expr, n, 'cos(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0) 
            # Check frequency
            if bb.is_constant() and abs(bb) > pi:
                print("Warning: Frequency may be out of range")    
                
            result = self.termXq(expr, n, k, q, lower, upper)
            # Make copies for the transformation:
            # X(q) = (exp(j*bb)*Xq (q*exp(j*bb)) + exp(-j*bb)*Xq(q*exp(-j*bb))) / 2 
            rq1 = deepcopy(result)
            rq2 = deepcopy(result)
            # Make general shifts 
            rq1.subs(q, q * sym.exp(sym.I * bb))           
            rq2.subs(q, q * sym.exp(-sym.I * bb))
            
            # Check special case shifts
            k0 = bb * self.N / 2 / pi
            if k0.is_integer and result.has_special:
                rq1.shift_k(k0)
                rq2.shift_k(-k0)
            else:
                rq1.rm_cases()
                rq2.rm_cases()                
            rq1.multiply(const * sym.exp(sym.I * cc) / 2)
            rq2.multiply(const * sym.exp(-sym.I * cc) / 2)
            
            # Add both parts 
            rq1.add(rq2)
            return rq1

        # Handle * exp(bb*n+cc)
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac):
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)                 
            result = self.termXq(expr, n, k, q, lower, upper)
            result.subs(q, q * sym.exp(bb))
            # No special cases remain
            result.rm_cases()
            result.multiply(const * sym.exp(cc))
            return result

        # Handle * a**n
        elif is_multiplied_with(expr, n, 'a**n', xn_fac):
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            lam = ref[0]
            bb = ref[1].coeff(n, 1)
            cc = ref[1].coeff(n, 0) 
            result = self.termXq(expr, n, k, q, lower, upper)
            result.subs(q, q * lam**bb)
            # No special cases remain
            result.rm_cases()
            result.multiply(const * lam**cc)           
            return result        

        # Handle * n       
        elif is_multiplied_with(expr, n, 'n', xn_fac):
            expr = expr / xn_fac[-1]
            result= self.termXq(expr, n, k, q, lower, upper)
            result.Xq = q * sym.diff(result.Xq, q) 
            if result.has_special:
                raise ValueError("No sym.diff possible for discrete cases, refine handles")
            else:
                return result                
        
        # No case matches so return None as this can be used to check if
        # a case was found and if a simplify is useful for the general
        # summation term simplify runs too long
        return None
    
    def term(self, expr, n, k):
       
        const, expr = factor_const(expr, n)
                
        if expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            result = const * self.function(expr, n, k)
            if self.is_inverse:
                result /= self.N
            return result  
        
        #
        # Transform of type X(k)
        #
        # Handle  b/(1 - a*exp(-j*2*pi*n/N))   
        # TODO can be extended to more complicated rational functions
        qq = sym.Symbol('qq')
        expr_qq = expr.replace(sym.exp(-sym.I * 2 * pi * n / self.N), qq)
        rat = expr_qq.as_numer_denom()
        if not rat[0].has(n) and not rat[0].has(qq) and not rat[1].has(n) and rat[1].has(qq) and rat[1].as_poly(qq).is_linear:
            bb = rat[0]
            a0 = sym.expand(rat[1]).coeff(qq, 1)
            c0 = sym.expand(rat[1]).coeff(qq, 0) 
            aa = -a0 / c0
            if self.is_inverse:
                # Do not remove this / 2, that helps simplification                
                result = bb * aa**k / (1 - aa**self.N) / c0 / 2
            else:
                result = bb * self.N * aa**(self.N - k) / (1 - aa**self.N) / c0 / 2 + bb * self.N * UnitImpulse(k) / c0 / 2
            return 2 * const * result 
            
        # Handle b/(1 - a*exp(+j*2*pi*n/N)) 
        expr_qq = expr.replace(sym.exp(sym.I * 2 * pi * n / self.N), qq)
        rat = expr_qq.as_numer_denom()
        if not rat[0].has(n) and not rat[0].has(qq) and not rat[1].has(n) and rat[1].has(qq) and rat[1].as_poly(qq).is_linear:
            bb = rat[0]
            a0 = sym.expand(rat[1]).coeff(qq, 1)
            c0 = sym.expand(rat[1]).coeff(qq, 0)
            aa = -a0 / c0
            if not self.is_inverse:
                result = bb * aa**k / (1 - aa**self.N) / c0 * self.N / 2
            else:
                result = bb * aa**(self.N - k) / (1 - aa**self.N) / c0 / 2 + bb * UnitImpulse(k) / c0 / 2           
            return 2 * const * result
        
        # 
        #          Transforms of type X(q)
        #
        q = sym.Symbol('q')
        k = -k if self.is_inverse else k
        
        # Call transform
        res = self.termXq(expr, n, k, q, 0, self.N - 1)
        # Find handles
        if not res is None:
            # Simplify q**N terms
            res.simp_qN(q, self.N)
            # Make final transform by putting all cases and terms together
            res.make_transform(self.N, self.is_inverse, q, sym.exp(-sym.I * 2 * pi / self.N * k), k)
            result = res.Xk
            
            # TODO smart simplification since simplify takes too long or
            # makes sometimes weird expressions might be better in
            # make_transform
            
            # result = sym.simplify(result)   #only useful if handles found
        # No handle found
        else:
            # Else makes it easy to exclude from simplifications 
            result = sym.summation(expr * sym.exp(-j * 2 * pi * n * k / self.N), (n, 0, self.N-1))
        
        result *= const
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


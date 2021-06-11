"""This module provides support for Laplace transforms.  It acts as a
quantity for SymPy's Laplace transform.  It calculates the unilateral
Laplace transform using:

F(s) = lim_{t_0\rightarrow 0} \int_{-t_0}^{\infty} f(t) e^{-s t} dt

In comparison, SymPy uses:

F(s) = \int_{0}^{\infty} f(t) e^{-s t} dt

The latter gives 0.5 for the Laplace transform of DiracDelta(t)
whereas the former version gives 1.  Note, SymPy is inconsistent in
that it gives DiracDelta(t) for the inverse Laplace transform of 1.

Another difference with this implementation is that it will transform
undefined functions such as v(t) to V(s).

These functions are for internal use by Lcapy.  

Copyright 2016--2021 Michael Hayes, UCECE

"""

from .transformer import UnilateralForwardTransformer
from .ratfun import Ratfun
from .sym import sympify, simplify, AppliedUndef
from .utils import factor_const, scale_shift, as_sum_terms
import sympy as sym

__all__ = ('LT', 'laplace_transform')

class LaplaceTransformer(UnilateralForwardTransformer):

    name = 'inverse Laplace transform'
    
    def noevaluate(self, expr, t, s):

        t0 = sympify('t0')
        return sym.Limit(sym.Integral(expr * sym.exp(-s * t),
                                      (t, t0, sym.oo)),
                         t0, 0, dir='-')

    def check(self, expr, t, s, **assumptions):

        if expr.has(s):
            self.error('Expression depends on s')

    def key(self, expr, t, s, **assumptions):
        return expr, t, s,

    def limits(self, expr, t, s, tmin, tmax):

        F = sym.integrate(expr * sym.exp(-s * t), (t, tmin, tmax))

        if not F.has(sym.Integral):
            return F

        if not F.is_Piecewise:
            self.error('Expecting piecewise')

        F, cond = F.args[0]
        if F.has(sym.Integral):
            self.error('Expecting integral')

        return F

    def do_0minus(self, expr, t, s):

        t0 = sym.symbols('t0', negative=True, real=True)

        F = self.limits(expr, t, s, t0, sym.oo)
        return sym.limit(F, t0, 0)

    def do_0(self, expr, t, s):

        return self.limits(expr, t, s, 0, sym.oo)

    def func(self, expr, t, s, inverse=False):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], t)    

        ssym = sympify(str(s))

        # Convert v(t) to V(s), etc.
        name = expr.func.__name__
        func = name[0].upper() + name[1:] + '(%s)' % s

        result = sympify(func).subs(ssym, s / scale) / abs(scale)

        if shift != 0:
            result = result * sym.exp(s * shift / scale)    
        return result


    def integral(self, expr, t, s):

        const, expr = factor_const(expr, t)

        if len(expr.args) != 2:
            self.error('Expecting two args')

        integrand = expr.args[0]

        if not isinstance(expr, sym.Integral):
            self.error('Expecting integral')

        if len(expr.args[1]) != 3:
            self.error('Require definite integral')

        var = expr.args[1][0]
        limits = expr.args[1][1:]
        const2, expr2 = factor_const(integrand, var)

        if (expr2.is_Function and
            expr2.args[0] == t - var and limits[0] == 0 and limits[1] == sym.oo):
            return const2 * self.term(expr2.subs(t - var, t), t, s) / s

        # Look for convolution integral
        # TODO, handle convolution with causal functions.
        if (limits[0] != -sym.oo) or (limits[1] != sym.oo):
            self.error('Need indefinite limits')

        if ((len(expr.args) != 2) or not expr2.is_Mul or
            not expr2.args[0].is_Function or not expr2.args[1].is_Function):
            self.error('Need integral of product of two functions')

        f1 = expr2.args[0]
        f2 = expr2.args[1]    
        # TODO: apply similarity theorem if have f(a * tau) etc.

        if (f1.args[0] == var and f2.args[0] == t - var):
            F1 = self.term(f1, var, s)
            F2 = self.term(f2.subs(t - var, t), t, s)
        elif (f2.args[0] == var and f1.args[0] == t - var):
            F1 = self.term(f1.subs(t - var, t), t, s)
            F2 = self.term(f2, var, s)
        else:            
            self.error('Cannot recognise convolution')

        return const2 * F1 * F2

    def derivative_undef(self, expr, t, s):

        if not isinstance(expr, sym.Derivative):
            self.error('Expecting derivative')

        if (not isinstance(expr.args[0], AppliedUndef) and
            expr.args[1][0] != t):
            self.error('Expecting function of t')

        ssym = sympify(str(s))    
        name = expr.args[0].func.__name__    
        func1 = name[0].upper() + name[1:] + '(%s)' % str(ssym)    
        return sympify(func1).subs(ssym, s) * s ** expr.args[1][1]

    def sin_cos(self, expr, t, s):

        # Sympy sometimes has problems with this...

        if not expr.is_Function or expr.func not in (sym.sin, sym.cos):
            self.error('Expression not sin or cos')
        arg = expr.args[0]
        a, b = scale_shift(arg, t)

        d = s**2 + a**2

        if expr.func is sym.sin:
            return (a * sym.cos(b) + s * sym.sin(b)) / d
        else:
            return (-a * sym.sin(b) + s * sym.cos(b)) / d        

    def term(self, expr, t, s):

        # Unilateral LT ignores expr for t < 0 so remove Piecewise.
        if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
            expr = expr.args[0].args[0]

        const, expr = factor_const(expr, t)

        terms = expr.as_ordered_terms()
        if len(terms) > 1:
            result = 0
            for term in terms:
                result += self.term(term, t, s)
            return const * result

        tsym = sympify(str(t))
        expr = expr.replace(tsym, t)

        if expr.has(sym.Integral):
            return self.integral(expr, t, s) * const

        if expr.is_Function and expr.func in (sym.sin, sym.cos) and expr.args[0].has(t):
            return self.sin_cos(expr, t, s) * const

        if expr.has(AppliedUndef):

            if expr.has(sym.Derivative):
                return self.derivative_undef(expr, t, s) * const    

            factors = expr.as_ordered_factors()
            if len(factors) == 1:
                return const * self.func(factors[0], t, s)
            elif len(factors) > 2:
                self.error('Cannot handle product')

            foo = factors[1]
            if foo.is_Function and foo.func == sym.exp and foo.args[0].has(t):
                scale, shift = scale_shift(foo.args[0], t)
                if shift == 0: 
                    result = self.func(factors[0], t, s)
                    return const * result.subs(s, s - scale)
            self.error('Cannot handle product')

        if expr.has(sym.Heaviside(t)):
            return self.do_0(expr.replace(sym.Heaviside(t), 1), t, s) * const

        if expr.has(sym.DiracDelta) or expr.has(sym.Heaviside):
            try:
                return self.do_0minus(expr, t, s) * const
            except ValueError:
                pass

        return self.do_0(expr, t, s) * const


laplace_transformer = LaplaceTransformer()


def laplace_transform(expr, t, s, evaluate=True, **assumptions):
    """Compute unilateral Laplace transform of expr with lower limit 0-.    
    """
    
    return laplace_transformer.transform(expr, t, s,
                                         evaluate=evaluate,
                                         **assumptions)


def LT(expr, t, s, evaluate=True, **assumptions):
    """Compute unilateral Laplace transform of expr with lower limit 0-.    
    """    

    return laplace_transformer.transform(expr, t, s,
                                         evaluate=evaluate,
                                         **assumptions)


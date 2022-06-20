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

Copyright 2016--2022 Michael Hayes, UCECE

"""

from .transformer import UnilateralForwardTransformer
from .ratfun import Ratfun
from .sym import sympify, simplify, AppliedUndef
from .utils import (factor_const, scale_shift, as_sum_terms, similarity_shift,
                    expand_functions)
from .extrafunctions import rect, tri, ramp, rampstep
import sympy as sym
from warnings import warn

__all__ = ('LT', 'laplace_transform')


class LaplaceTransformer(UnilateralForwardTransformer):

    name = 'Laplace transform'

    def noevaluate(self, expr, t, s):

        t0 = sympify('t0')
        return sym.Limit(sym.Integral(expr * sym.exp(-s * t),
                                      (t, t0, sym.oo)),
                         t0, 0, dir='-')

    def check(self, expr, t, s, **kwargs):

        if expr.has(s):
            self.error('Expression depends on s')

    def key(self, expr, t, s, **kwargs):
        return (expr, t, s,
                kwargs.get('zero_initial_conditions', True))
        return expr, t, s,

    def integrate(self, expr, t, s, tmin, tmax):

        F = sym.integrate(expr * sym.exp(-s * t), (t, tmin, tmax))

        if not F.has(sym.Integral):
            return F

        if isinstance(F, sym.Integral):
            self.error('Integration failed')

        if not F.is_Piecewise:
            self.error('Expecting piecewise')

        F, cond = F.args[0]
        if F.has(sym.Integral):
            self.error('Expecting integral')

        return F

    def integrate_0minus(self, expr, t, s):

        t0 = sym.symbols('t0', negative=True, real=True)

        F = self.integrate(expr, t, s, t0, sym.oo)
        return sym.limit(F, t0, 0)

    def integrate_0(self, expr, t, s):

        return self.integrate(expr, t, s, 0, sym.oo)

    def func(self, expr, t, s, inverse=False):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], t)

        # Convert v(t) to V(s), etc.
        name = expr.func.__name__
        func = sym.Function(name[0].upper() + name[1:])

        result = func(s / scale) / abs(scale)

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

        # Handle integral(x(tau), (tau, -oo, t))
        if not limits[0].is_positive and limits[1] == t:
            if isinstance(expr2, AppliedUndef):
                if expr2.args[0] == expr.args[1][0]:
                    return const2 * self.func(expr2, expr2.args[0], s) / s

        # Look for convolution integral
        if limits[0].is_positive:
            self.error('Cannot handle lower limit %s' % limits[0])

        if limits[1] < t:
            self.error('Cannot handle upper limit %s' % limits[1])

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

    def derivative_undef(self, expr, t, s, zero_initial_conditions=True):

        if not isinstance(expr, sym.Derivative):
            self.error('Expecting derivative')

        if (not isinstance(expr.args[0], AppliedUndef) and
                expr.args[1][0] != t):
            self.error('Expecting function of t')

        name = expr.args[0].func.__name__
        func1 = sym.Function(name[0].upper() + name[1:])
        order = expr.args[1][1]
        result = func1(s) * s ** order

        if not zero_initial_conditions:
            # Handle initial conditions.  FIXME  use 0^- for 0.
            v = sym.Function(name)(t)
            for m in range(order):
                result -= s**(order - m - 1) * \
                    sym.Derivative(v, t, m).subs(t, 0)
        return result

    def sin_cos(self, expr, t, s):

        # Handle exp(-alpha * t) * sin(omega * t + phi) * u(t - tau)
        # The exp(-alpha * t) and u(t - tau) parts are optional.

        # Sympy sometimes has problems with this...

        factors = expr.as_ordered_factors()

        if len(factors) > 3:
            raise ValueError('Not expsin, too many factors')

        alpha = 0
        beta = 0
        m = 0
        if (factors[m].is_Function and factors[m].func is sym.exp):
            exparg = factors[m].args[0]
            alpha, beta = scale_shift(exparg, t)
            m += 1

        if not (factors[m].is_Function and factors[m].func in (sym.sin, sym.cos)):
            raise ValueError('Not expsin, no sin/cos')

        sincosarg = factors[m].args[0]
        omega, phi = scale_shift(sincosarg, t)

        if factors[m].func is sym.cos:
            phi += sym.pi / 2

        m += 1

        if len(factors) == m + 1 and not (factors[m].is_Function and factors[m].func is sym.Heaviside):
            raise ValueError('Not expsin, no Heaviside')

        tau = 0
        if len(factors) == m + 1:
            eta, zeta = scale_shift(factors[m].args[0], t)
            if eta != 1:
                raise ValueError('Need to use similarity theorem')
            tau = -zeta

            if tau.is_negative:
                tau = 0
            elif tau.is_Symbol and not tau.is_positive:
                print('Assuming %s is positive' % tau)
            elif (tau.is_Mul and len(tau.args) == 2 and
                  tau.args[0].is_negative and tau.args[1].is_Symbol):
                print('Assuming %s is positive' % -tau)
                tau = 0

        if tau != 0:
            phi += omega * tau

        E = (omega * sym.cos(phi) + (s - alpha) *
             sym.sin(phi)) / (omega**2 + (s - alpha)**2)

        if tau != 0:
            E *= sym.exp(-tau * s)
            if alpha != 0:
                E *= sym.exp(alpha * tau)

        if beta != 0:
            E = sym.exp(beta) * E

        return E

    def function(self, expr, t, s):

        expr, scale, shift = similarity_shift(expr, t)

        if shift != 0:
            # This will be handled by rewriting the function.
            return None

        if expr.func is rect:
            warn('Laplace transform ignores rect(t) for t < 0')
            return (1 - sym.exp(-s / (2 * scale))) / s

        elif expr.func is tri:
            warn('Laplace transform ignores tri(t) for t < 0')
            return 1 / s - (1 - sym.exp(-s / scale)) / (scale * s**2)

        elif expr.func is ramp:
            return scale / s**2

        elif expr.func is rampstep:
            return (1 - sym.exp(-s / scale)) / s**2

        return None

    def term(self, expr, t, s, **kwargs):

        # Unilateral LT ignores expr for t < 0 so remove Piecewise.
        if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
            expr = expr.args[0].args[0]

        if expr.has(sym.sinh, sym.cosh, sym.tanh):
            expr = expr.rewrite(sym.exp)

        const, expr = factor_const(expr, t)

        terms = expr.expand(deep=False).as_ordered_terms()
        if len(terms) > 1:
            result = 0
            for term in terms:
                result += self.term(term, t, s, **kwargs)
            return const * result

        tsym = sympify(str(t))
        expr = expr.replace(tsym, t)

        if expr.has(sym.Integral):
            return self.integral(expr, t, s) * const

        if expr.has(sym.sin, sym.cos):
            try:
                return self.sin_cos(expr, t, s) * const
            except:
                pass

        if expr.has(AppliedUndef):

            if (expr.is_Mul and len(expr.args) == 2 and
                isinstance(expr.args[0], sym.DiracDelta) and
                    isinstance(expr.args[1], (AppliedUndef, sym.Subs))):
                return expr.args[1]

            if expr.has(sym.Derivative):
                return self.derivative_undef(expr, t, s, **kwargs) * const

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

        if expr.is_Function and expr.args[0].has(t):
            result = self.function(expr, t, s)
            if result is not None:
                return result * const
            expr = expand_functions(expr, t)

        if expr.has(sym.Heaviside(t)):
            return self.integrate_0(expr.replace(sym.Heaviside(t), 1), t, s) * const

        if expr.has(sym.DiracDelta) or expr.has(sym.Heaviside):
            try:
                return self.integrate_0minus(expr, t, s) * const
            except ValueError:
                pass

        return self.integrate_0(expr, t, s) * const


laplace_transformer = LaplaceTransformer()


def laplace_transform(expr, t, s, evaluate=True, **kwargs):
    """Compute unilateral Laplace transform of expr with lower limit 0-.
    """

    return laplace_transformer.transform(expr, t, s,
                                         evaluate=evaluate,
                                         **kwargs)


def LT(expr, t, s, evaluate=True, **kwargs):
    """Compute unilateral Laplace transform of expr with lower limit 0-.
    """

    return laplace_transformer.transform(expr, t, s,
                                         evaluate=evaluate,
                                         **kwargs)

"""This module provides support for Fourier transforms.  It calculates
the bilateral Fourier transform using:

   S(f) = \int_{-\infty}^{\infty} s(t) e^{-j * 2 * \pi * t} dt

It also allows functions that strictly do not have a Fourier transform
by using Dirac deltas.  For example, a, cos(a * t), sin(a * t), exp(j
* a * t).


Copyright 2016--2022 Michael Hayes, UCECE

"""

# TODO:
# Add DiracDelta(t, n)
# Simplify  (-j * DiracDelta(f - 1) + j * DiracDelta(f + 1)).inverse_fourier()
# This should give 2 * sin(2 * pi * t)

from sympy.core.function import AppliedUndef
from sympy import sympify, pi, exp, I, oo, S, sign, sin, cos, sinh, cosh, tanh
from sympy import DiracDelta, Heaviside, FourierTransform, Integral
from sympy import fourier_transform as sympy_fourier_transform, Function
from .sym import symsimplify, j
from .transformer import BilateralForwardTransformer
from .utils import factor_const, similarity_shift, expand_functions
from .extrafunctions import rect, sincn, sincu, trap, tri

__all__ = ('FT', 'IFT')


class FourierTransformer(BilateralForwardTransformer):

    name = 'Fourier transform'

    def key(self, expr, t, f, **kwargs):
        return expr, t, f

    def simplify_term(self, expr, var):

        return symsimplify(expr)

    def sympy(self, expr, t, f):

        result = sympy_fourier_transform(expr, t, f)
        if expr != 0 and result == 0:
            # There is a bug in SymPy where it returns 0.
            self.error()

        if isinstance(result, FourierTransform):
            self.error(
                'Try using partfrac() for partial fraction expansion first.')

        return result

    def func(self, expr, t, f):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        # Convert v(t) to V(f), etc.
        name = expr.func.__name__
        if self.is_inverse:
            func = Function(name[0].lower() + name[1:])
        else:
            func = Function(name[0].upper() + name[1:])

        result = func(f)
        return result

    def integral(self, expr, t, f):

        const, expr = factor_const(expr, t)

        if len(expr.args) != 2:
            self.error()

        integrand = expr.args[0]

        if not isinstance(expr, Integral):
            self.error()

        if len(expr.args[1]) != 3:
            self.error('Require definite integral')

        var = expr.args[1][0]
        limits = expr.args[1][1:]
        const2, expr2 = factor_const(integrand, var)

        if (expr2.is_Function and
                expr2.args[0] == t - var and limits[0] == 0 and limits[1] == oo):
            return const2 * self.term(expr2.subs(t - var, t), t, f) / f

        # Look for convolution integral
        # TODO, handle convolution with causal functions.
        if (limits[0] != -oo) or (limits[1] != oo):
            self.error('Need indefinite limits')

        if ((len(expr.args) != 2) or not expr2.is_Mul or
                not expr2.args[0].is_Function or not expr2.args[1].is_Function):
            self.error('Need integral of product of two functions')

        f1 = expr2.args[0]
        f2 = expr2.args[1]
        # TODO: apply similarity theorem if have f(a * tau) etc.

        if (f1.args[0] == var and f2.args[0] == t - var):
            F1 = self.term(f1, var, f)
            F2 = self.term(f2.subs(t - var, t), t, f)
        elif (f2.args[0] == var and f1.args[0] == t - var):
            F1 = self.term(f1.subs(t - var, t), t, f)
            F2 = self.term(f2, var, f)
        else:
            self.error('Cannot recognise convolution')

        return const2 * F1 * F2

    def function(self, expr, t, f):

        # Handle expressions with a function of FOO, e.g.,
        # v(t), v(t) * y(t),  3 * v(t) / t, v(4 * a * t), etc.,

        if not expr.has(AppliedUndef):
            self.error()

        const, expr = factor_const(expr, t)

        if isinstance(expr, AppliedUndef) and expr.args[0] == t:
            return self.func(expr, t, f) * const

        tsym = sympify(str(t))
        expr = expr.subs(tsym, t)

        rest = S.One
        undefs = []
        for factor in expr.as_ordered_factors():
            if isinstance(factor, AppliedUndef):
                if factor.args[0] != t:
                    self.error('Weird function %s not of %s' % (factor, t))
                undefs.append(factor)
            else:
                rest *= factor

        if rest.has(AppliedUndef):
            # Have something like 1/v(t)
            self.error()

        exprs = undefs
        if rest.has(t):
            exprs = exprs + [rest]
            rest = S.One

        result = self.term(exprs[0], t, f) * rest

        if len(exprs) == 1:
            return result * const

        for m in range(len(exprs) - 1):
            nu = self.dummy_var(expr, 'tau' if self.is_inverse else 'nu',
                                level=m, real=True)
            expr2 = self.term(exprs[m + 1], t, f)
            result = Integral(result.subs(f, f - nu) * expr2.subs(f, nu),
                              (nu, -oo, oo))

        return result * const

    def term(self, expr, t, f):

        const, expr = factor_const(expr, t)

        if expr.has(Integral):
            return self.integral(expr, t, f) * const

        if isinstance(expr, AppliedUndef) and expr.args[0] == t:
            return self.func(expr, t, f) * const

        # TODO add u(t) <-->  delta(f) / 2 - j / (2 * pi * f)

        if expr.has(AppliedUndef) and expr.args[0] == t:
            # Handle v(t), v(t) * y(t),  3 * v(t) / t etc.
            return self.function(expr, t, f) * const

        # Check for constant.
        if not expr.has(t):
            return expr * DiracDelta(f) * const

        one = S.One
        const1 = const
        other = one
        exps = one
        factors = expr.expand().as_ordered_factors()
        for factor in factors:
            if not factor.has(t):
                const1 *= factor
            else:
                if factor.is_Function and factor.func == exp and factor.args[0].has(I):
                    exps *= factor
                else:
                    other *= factor

        sf = -f if self.is_inverse else f

        if other != 1 and exps == 1:
            if other == t:
                return const1 * I / (2 * pi) * DiracDelta(sf, 1)
            elif other == t**2:
                return -const1 / (2 * pi)**2 * DiracDelta(sf, 2)
            # TODO check for other powers of t...
            elif other == sign(t):
                return const1 / (I * pi * sf)
            elif other == sign(t) * t:
                return -const1 * 2 / (2 * pi * f)**2
            elif other == Heaviside(t):
                return const1 / (I * 2 * pi * f) + const1 * DiracDelta(sf) / 2
            elif other == 1 / t:
                return -const1 * I * pi * sign(sf)
            elif other == 1 / t**2:
                return -const1 * 2 * pi**2 * sf * sign(sf)
            elif False and other == exp(t):
                # SymPy does not like Dirac Delta of a complex argument.
                return DiracDelta(sf + I / (2 * pi))
            elif other.is_Function and other.func == Heaviside and other.args[0] == t:
                return (const1 / (I * 2 * pi * sf) + const1 * DiracDelta(sf) / 2)
            elif other == Heaviside(t) * t:
                return -const1 / (2 * pi * f)**2 + const1 * I * DiracDelta(sf, 1) / (4 * pi)
            # t * u(t - tau)
            elif (other.is_Mul and len(other.args) == 2 and
                  other.args[0] == t and other.args[1].is_Function and
                  other.args[1].func == Heaviside and other.args[1].args[0] == t):
                e = exp(I * 2 * pi * sf)
                return I * DiracDelta(sf, 1) / (4 * pi) * e - 1 / (4 * pi**2 * f**2) * e
            # exp(c1 * t + c0) * u(t)
            elif (other.is_Mul and len(other.args) == 2 and
                  other.args[1].is_Function and other.args[1].func == exp and
                  other.args[1].args[0].is_polynomial(t) and
                  other.args[1].args[0].as_poly(t).is_linear and
                  other.args[0].is_Function and
                  other.args[0].func == Heaviside and
                  other.args[0].args[0] == t):
                foo = other.args[1].args[0]
                c0 = foo.coeff(t, 0)
                c1 = foo.coeff(t, 1)
                return const1 * exp(c0) / (I * 2 * pi * sf - c1)
            elif other.is_Function and other.func == sincn and other.args[0] == t:
                return const1 * rect(f)
            elif other.is_Function and other.func == sincu and other.args[0] == t:
                return const1 * pi * rect(f * pi)
            elif (other.is_Pow and other.args[1] == 2 and
                  other.args[0].is_Function and other.args[0].func == sincn and
                  other.args[0].args[0] == t):
                other = other.args[0]
                return const1 * tri(f)
            elif other.is_Function and other.func == rect and other.args[0] == t:
                return const1 * sincn(f)
            elif other.is_Function and other.func == tri and other.args[0] == t:
                return const1 * sincn(f)**2
            elif other.is_Function and other.func == trap and other.args[0] == t:
                alpha = other.args[1]

                # Check for rect
                if alpha == 0:
                    return const1 * sincn(f)
                return alpha * const1 * sincn(f) * sincn(alpha * f)

            #
            # factor = other.factor()
            # const2, factor2 = factor_const(factor, t)
            # if factor2.is_Pow and factor2.args[1] == -2:
            #     foo = factor2.args[0]
            #     a = foo.coeff(t, 1)
            #     b = foo.coeff(t, 0)
            #     if a != 0:
            #         return const1 * const2 * f * exp(-b * 2 *pi * f / a) * sign(f) / (2 * pi)

            # Sympy incorrectly gives exp(-a * t) instead of exp(-a * t) * Heaviside(t)
            elif other.is_Pow and other.args[1] == -1 and other.args[0].has(t):
                foo = other.args[0]

                if (foo.is_polynomial(t) and foo.as_poly(t).is_linear and
                        foo.is_complex):
                    c0 = foo.coeff(t, 0)
                    c1 = foo.coeff(t, 1)
                    s = (2 * pi * I) / c1
                    return const1 * s * exp(c0 * sf * s) * Heaviside(-sf)
                elif foo.is_Function and foo.func == cosh and foo.args[0] == t:
                    return const * pi / cosh(pi**2 * sf)
                elif foo.is_Function and foo.func == sinh and foo.args[0] == t:
                    return -I * const * pi * tanh(pi**2 * sf)
            elif other.is_Function and other.func == tanh and other.args[0] == t:
                return -I * const * pi / sinh(pi**2 * sf)

            if expr == t * DiracDelta(t, 1):
                return const * sf / (-I * 2 * pi)

            # Apply similarity and shift theorems.
            expr2, scale, shift = similarity_shift(expr, t)
            if scale != 1 or shift != 0:
                result = self.term(expr2, t, f / scale) / abs(scale)
                if shift != 0:
                    result *= exp(-I * 2 * pi * f / scale * shift)
                return const * result

            if expr.is_Function and expr.args[0] == t:
                # Try to handle functions such as ramp, rampstep.
                expr2 = expand_functions(expr, t)
                if expr != expr2:
                    terms = expr2.as_ordered_terms()
                    result = 0
                    for term in terms:
                        result += self.term(term, t, f)
                    return result * const

            # Punt and use SymPy.  Should check for t**n, t**n * exp(-a * t), etc.
            return const * self.sympy(expr, t, sf)

        args = exps.args[0]
        foo = args / t
        if foo.has(t):
            # Have exp(a * t**n), SymPy might be able to handle this
            return const * self.sympy(expr, t, sf)

        if exps != 1 and foo.has(I):
            if other == 1:
                return const1 * DiracDelta(sf - foo / (I * 2 * pi))
            Q = self.term(other, t, f)
            return const1 * Q.subs(f, (f - foo / (I * 2 * pi)))

        return const * self.sympy(expr, t, sf)

    def noevaluate(self, expr, t, f):

        return Integral(expr * exp(-j * 2 * pi * f * t), (t, -oo, oo))

    def check(self, expr, t, f):

        if expr.has(f):
            self.error('Expression depends on f')

        if expr.is_Piecewise and expr.args[0].args[1].has(t >= 0):
            self.error('Expression is unknown for t < 0 (use causal=True)')

    def rewrite(self, expr, var):

        # sym.rewrite(exp) can create exp(log...)
        if expr.has(sin):
            expr = expr.replace(lambda expr: expr.is_Function and expr.func == sin,
                                lambda expr: expr.rewrite(exp))
        if expr.has(cos):
            expr = expr.replace(lambda expr: expr.is_Function and expr.func == cos,
                                lambda expr: expr.rewrite(exp))

        expr = expr.expand()
        return expr


fourier_transformer = FourierTransformer()


def FT(expr, t, f, **kwargs):
    """Compute bilateral Fourier transform of expr.

    Undefined functions such as v(t) are converted to V(f)

    This also handles some expressions that do not really have a Fourier
    transform, such as a, cos(a * t), sin(a * t), exp(I * a * t)."""

    return fourier_transformer.transform(expr, t, f, **kwargs)


def fourier_transform(expr, t, f, **kwargs):
    """Compute bilateral Fourier transform of expr.

    Undefined functions such as v(t) are converted to V(f)

    This also handles some expressions that do not really have a Fourier
    transform, such as a, cos(a * t), sin(a * t), exp(I * a * t)."""

    return fourier_transformer.transform(expr, t, f, **kwargs)

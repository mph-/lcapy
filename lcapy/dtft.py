"""This module provides support for discrete-time Fourier transforms.

 It calculates the discrete-time Fourier transform using:

   X(f) = \sum_{n=\infty}^{\infty} x(n) e^{-j * 2 * \pi * n * dt * f}

Copyright 2021--2022 Michael Hayes, UCECE

"""

import sympy as sym
from sympy import oo, DiracDelta
from sympy.core import S
from .transformer import BilateralForwardTransformer
from .sym import sympify, AppliedUndef, j, pi, miscsymbol
from .sym import dt
from .utils import factor_const, scale_shift
from .ztransform import is_multiplied_with
from .extrafunctions import UnitImpulse, UnitStep, sincu, sincn, dtrect, dtsign, tri
from warnings import warn


__all__ = ('DTFT', )


class DTFTTransformer(BilateralForwardTransformer):

    name = 'DTFT'

    def key(self, expr, n, f, **kwargs):
        return expr, n, f, kwargs.get('images', 0)

    def noevaluate(self, expr, n, f):

        foo = expr * sym.exp(-2 * j * pi * n * dt * f)
        result = sym.Sum(foo, (n, -oo, oo))
        return result

    def check(self, expr, n, f, images=0, **kwargs):

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
            self.error('Expression is unknown for n < 0 (use causal=True)')

    def add_images(self, expr, f):
        if self.m1 == self.m2:
            return expr

        msym = self.dummy_var(expr, 'm', level=0, integer=True)
        foo = expr.replace(f, f - msym / dt)
        result = sym.Sum(foo, (msym, self.m1, self.m2))
        return result

    def sympy(self, expr, n, f):

        foo = expr * sym.exp(-2 * j * pi * n * dt * f)
        result = sym.summation(foo, (n, -oo, oo))

        return result

    def func(self, expr, n, f):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], n)

        # Convert v(n) to V(f), etc.
        name = expr.func.__name__
        func = sym.Function(name[0].upper() + name[1:])

        result = func(f / scale) / abs(scale)

        if shift != 0:
            result = result * sym.exp(2 * sym.I * sym.pi * f * shift / scale)
        # Perhaps return X_(1/dt)(f) but how to denote?
        return self.add_images(result, f)

    def function(self, expr, n, f, **kwargs):

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

        result = self.term(exprs[0], n, f, **kwargs) * rest

        if len(exprs) == 1:
            return result * const

        self.error('TODO')

    def term(self, expr, n, f, **kwargs):

        const, expr = factor_const(expr, n)
        args = expr.args
        twopidt = 2 * sym.pi * dt
        xn_fac = []

        # Check for constant.
        if not expr.has(n):
            return self.add_images(expr * DiracDelta(f) * const, f) / dt

        elif expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            return self.function(expr, n, f, **kwargs) * const

        # Handle step u(n-n0)  or u(-n-n0)  only
        elif expr.is_Function and expr.func in (sym.Heaviside, UnitStep):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)
            delay = -bb / aa
            if delay.is_integer or delay.is_integer is None:
                res = 1 / (1 - sym.exp(-sym.I * twopidt * f))
                if aa.is_negative:
                    res = res.subs(f, -f)
                return const * sym.exp(-sym.I * twopidt * f * delay) * res + const * self.add_images(DiracDelta(f), f) / dt / 2

        # Handle exp(j*a*n+b) * x(n)    o--o   X(W-a) * exp(b)
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac) and abs(xn_fac[-1] / sym.exp(args[0].coeff(n, 0))) == 1:
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            aa = ref[0].coeff(n, 1) / sym.I
            bb = ref[0].coeff(n, 0)
            X = self.term(expr, n, f, **kwargs)
            res = X.subs(f,  f - aa / twopidt)
            return const * res * sym.exp(bb)

        # Handle sin(b*n+c) * x(n)    o--o   j/2 (exp(-jc) * X(W+b) - X(W-b) * exp(jc))
        elif is_multiplied_with(expr, n, 'sin(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)
            X = self.term(expr, n, f, **kwargs)
            Xp = X.subs(f, f + bb / twopidt)
            Xm = X.subs(f, f - bb / twopidt)
            res = sym.I * (Xp * sym.exp(-sym.I * cc) -
                           Xm * sym.exp(sym.I * cc))
            return const / 2 * res

        # Handle cos(b*n+c)*x(n)    o--o   1/2 (exp(-jc)* X(W+b) + X(W-b) * exp(jc))
        elif is_multiplied_with(expr, n, 'cos(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)
            X = self.term(expr, n, f, **kwargs)
            Xp = X.subs(f, f + bb / twopidt)
            Xm = X.subs(f, f - bb / twopidt)
            res = (Xp * sym.exp(-sym.I * cc) + Xm * sym.exp(sym.I * cc))
            return const / 2 * res

        # Multiplication with n       use n * x(n)  o--o  j / twopidt * d/df X(f)
        elif is_multiplied_with(expr, n, 'n', xn_fac):
            expr = expr / xn_fac[-1]
            X = self.term(expr, n, f, **kwargs)
            return const / twopidt * sym.I * sym.simplify(sym.diff(X, f))

        # Handle u(n+n0) * a **n * exp(b*n+c)
        elif is_multiplied_with(expr, n, 'UnitStep', xn_fac):
            #
            expr /= xn_fac[-1]
            # Handle aa*n + bb  as argument of step
            ref = xn_fac[-1].args
            aa = ref[0].coeff(n, 1)
            bb = ref[0].coeff(n, 0)
            # Check argument of step function
            if abs(aa) != 1 or not (bb.is_integer or bb.is_integer is None):
                warn("Check argument of Step function")
            delay = -bb / aa
            #
            prefac = sym.exp(-sym.I * twopidt * f * delay)
            # Rename  summation index  aa*n+bbk = k --> n
            expr = expr.subs(n, aa * n + delay)

            # Handle u(n) * a **n * exp(b*n+c)
            e_fac = 1
            # Collect all exp. factors
            while is_multiplied_with(expr, n, 'exp(n)', xn_fac):
                expr /= xn_fac[-1]
                expr = sym.simplify(expr)
                ref = xn_fac[-1].args
                bb = ref[0].coeff(n, 1)
                cc = ref[0].coeff(n, 0)
                e_fac *= sym.exp(bb)
                prefac *= sym.exp(cc)

            # Collect all a factors
            while is_multiplied_with(expr, n, 'a**n', xn_fac):
                expr /= xn_fac[-1]
                expr = sym.simplify(expr)
                ref = xn_fac[-1].args
                lam = ref[0]
                bb = ref[1].coeff(n, 1)
                cc = ref[1].coeff(n, 0)
                e_fac *= lam ** bb
                prefac *= lam ** cc

            if not expr.has(n):
                if e_fac.is_integer and abs(e_fac) > 1:
                    warn("Check for convergence")
                res = expr / (1 - e_fac * sym.exp(-sym.I * twopidt * f))
                if aa.is_negative:
                    res = res.subs(f, -f)
                return const * prefac * res

        # Handle impulse delta (n-n0)
        elif expr.is_Function and expr.func == UnitImpulse:
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)
            delay = -bb / aa
            if delay.is_integer or delay.is_integer is None:
                return const * sym.exp(-sym.I * delay * twopidt * f)

        # Handle signum
        elif (len(args) == 1 and expr.is_Function and
              expr.func == dtsign and args[0].is_polynomial(n) and args[0].as_poly(n).is_linear):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)
            delay = -bb / aa
            if delay.is_integer:
                res = 1 / (1 - sym.exp(-sym.I * twopidt * f))
                if aa.is_negative:
                    res = res.subs(f, -f)
                return 2 * const * sym.exp(-sym.I * twopidt * f * delay) * res

        # Handle sincu/sincn
        elif (len(args) == 1 and expr.is_Function and
              (expr.func == sincu or expr.func == sincn) and
              args[0].is_polynomial(n) and args[0].as_poly(n).is_linear):
            aa = args[0].coeff(n, 1)
            if aa.is_number and abs(aa) > pi:
                warn("Argument out of range (-pi, pi)")
            bb = args[0].coeff(n, 0)
            delay = bb / aa
            K = aa / dt
            prefac = 1 / aa
            if expr.func == sincu:
                prefac *= sym.pi
                K /= sym.pi
            if delay.is_integer:
                result = const * prefac * \
                    dtrect(f / K) * sym.exp(sym.I * delay * twopidt * f)
                return self.add_images(result, f)

        # Handle sincu**2
        elif (expr.is_Pow and args[1] == 2 and
              args[0].is_Function and args[0].func in (sincu, sincn) and
              args[0].args[0].is_polynomial(n) and args[0].args[0].as_poly(n).is_linear):
            aa = args[0].args[0].coeff(n, 1)
            if aa.is_number and abs(aa) > pi:
                warn("Argument out of range (-pi, pi)")
            bb = args[0].args[0].coeff(n, 0)
            delay = bb / aa
            K = aa / dt
            prefac = 1 / aa
            if args[0].func == sincu:
                prefac *= sym.pi
                K *= sym.pi
            if delay.is_integer:
                result = const * prefac * \
                    tri(f / K) * sym.exp(sym.I * delay * twopidt * f)
                return self.add_images(result, f)

        # Handle dtrect
        elif (len(args) == 1 and expr.is_Function and expr.func == dtrect and
              args[0].is_polynomial(n) and args[0].as_poly(n).is_linear):
            N = 1 / args[0].coeff(n, 1)
            if N.is_negative:
                warn("Negative N for dtrect((n-n0)/N)")
            else:
                delay = N * args[0].coeff(n, 0)
                if delay.is_integer:
                    if N.is_even:
                        delay += S.Half
                    elif not N.is_odd:
                        if N.is_symbol:
                            warn(
                                "Assuming %s odd; if even use %s = symbol('%s', even=True)" % (N, N, N))
                        else:
                            warn("Assuming %s odd" % N)

                    return const * sym.exp(sym.I * delay * twopidt * f) * sym.sin(twopidt * f * N / 2) / sym.sin(twopidt * f / 2)

        return const * self.sympy(expr, n, f)


dtft_transformer = DTFTTransformer()


def discrete_time_fourier_transform(expr, n, f, images=0, evaluate=True,
                                    **kwargs):
    """Compute bilateral discrete-time Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(f)
    """

    return dtft_transformer.transform(expr, n, f, evaluate=evaluate,
                                      images=images, **kwargs)


def DTFT(expr, n, f, images=0, evaluate=True, **kwargs):
    """Compute bilateral discrete-time Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(f)
    """

    return dtft_transformer.transform(expr, n, f, evaluate=evaluate,
                                      images=images, **kwargs)

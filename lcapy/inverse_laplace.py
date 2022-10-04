"""This module provides support for the unilateral inverse Laplace
transform.

Copyright 2021--2022 Michael Hayes, UCECE

"""

from .transformer import UnilateralInverseTransformer
from .ratfun import Ratfun
from .sym import simplify, AppliedUndef
from .utils import factor_const, scale_shift, as_sum_terms
from .matrix import matrix_inverse
from warnings import warn
import sympy as sym

Zero = sym.S.Zero
One = sym.S.One

__all__ = ('ILT', 'inverse_laplace_transform')


class InverseLaplaceTransformer(UnilateralInverseTransformer):

    name = 'inverse Laplace transform'

    def noevaluate(self, expr, s, t):

        # Construct Bromwich integral.
        self.error('TODO')

    def check(self, expr, s, t, **kwargs):

        if expr.has(t):
            for term in expr.as_ordered_terms():
                # Allow Subs(Derivativee(undef, t), t, 0)
                if term.has(t) and not term.is_constant:
                    self.error('Expression depends on t')

    def key(self, expr, s, t, **kwargs):
        return (expr, s, t,
                kwargs.get('causal', False),
                kwargs.get('zero_initial_conditions', True),
                kwargs.get('damped_sin', True),
                kwargs.get('damping', None))

    def func(self, expr, s, t):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], s)

        # Convert V(s) to v(t), etc.
        name = expr.func.__name__
        undef = sym.Function(name[0].lower() + name[1:])(t)
        result = undef.subs(t, t / scale) / abs(scale)

        if shift != 0:
            result = result * sym.exp(t * shift / scale)
        return result

    def do_damped_sin(self, expr, s, t):

        ncoeffs, dcoeffs = expr.coeffs()
        K = ncoeffs[0] / dcoeffs[0]

        ncoeffs = [(c / ncoeffs[0]) for c in ncoeffs]
        dcoeffs = [(c / dcoeffs[0]) for c in dcoeffs]

        if len(ncoeffs) > 3 or len(dcoeffs) > 3:
            self.error('Not a second-order response')

        omega0 = sym.sqrt(dcoeffs[2])
        zeta = dcoeffs[1] / (2 * omega0)

        if zeta.is_constant() and zeta > 1:
            warn('Expression is overdamped')

        sigma1 = (zeta * omega0).simplify()
        omega1 = (omega0 * sym.sqrt(1 - zeta**2)).simplify()
        K = (K / omega1).simplify()

        E = sym.exp(-sigma1 * t)
        S = sym.sin(omega1 * t)

        h = K * E * S

        # If overdamped
        #h = K * sym.exp(-sigma1 * t) * sym.sinh(omega0 * mu * t)

        if len(ncoeffs) == 1:
            return Zero, h

        C = sym.cos(omega1 * t)
        kCd = omega1
        kSd = -sigma1
        hd = K * E * (kCd * C + kSd * S)

        if len(ncoeffs) == 2:
            return Zero, K * E * (kCd * C + (ncoeffs[1] + kSd) * S)

        kCdd = -2 * omega1 * sigma1
        kSdd = sigma1**2 - omega1**2

        G = K * E * ((kCdd + ncoeffs[1] * kCd) * C +
                     (kSdd + ncoeffs[1] * kSd + ncoeffs[2]) * S)

        return K * kCd * sym.DiracDelta(t), G

    def ratfun(self, expr, s, t, **kwargs):

        if kwargs.pop('pdb', False):
            import pdb
            pdb.set_trace()

        sexpr = Ratfun(expr, s)

        if kwargs.get('damped_sin', False):
            if sexpr.degree == 2:
                return self.do_damped_sin(sexpr, s, t)
            # if False and sexpr.degree == 3 and Ratfun(expr * s).degree == 2:
            #    return self.do_damped_sin3(sexpr, s, t)

        self.debug('Finding QRPO representation')

        damping = kwargs.get('damping', None)
        Q, R, P, O, delay, undef = sexpr.as_QRPO(damping)

        if delay != 0:
            # This will be caught and trigger expansion of the expression.
            self.error('Unhandled delay %s' % delay)

        cresult = Zero

        if Q:
            Qpoly = sym.Poly(Q, s)
            C = Qpoly.all_coeffs()
            for n, c in enumerate(C):
                cresult += c * sym.diff(sym.DiracDelta(t), t, len(C) - n - 1)

        if R == []:
            return cresult, 0

        uresult = 0
        for m, (p, n, A) in enumerate(zip(P, O, R)):

            # This is zero for the conjugate pole.
            if R[m] == 0:
                continue

            # Search and remove conjugate pair.
            has_conjpair = False
            if p.is_complex and kwargs.get('pairs', True):
                pc = p.conjugate()
                Ac = A.conjugate()
                for m2, p2 in enumerate(P[m + 1:]):
                    m2 += m + 1
                    if n == O[m2] and p2 == pc and R[m2] == Ac:
                        R[m2] = 0
                        has_conjpair = True
                        break

            if has_conjpair:
                # Combine conjugate pairs.
                p = p.expand(complex=True)
                A = A.expand(complex=True)
                p_re = sym.re(p)
                p_im = sym.im(p)
                A_re = sym.re(A)
                A_im = sym.im(A)
                et = sym.exp(p_re * t)
                result = 2 * A_re * et * sym.cos(p_im * t)
                result -= 2 * A_im * et * sym.sin(p_im * t)
            else:
                result = A * sym.exp(p * t)

            if n > 1:
                result *= t ** (n - 1) / sym.factorial(n - 1)
            uresult += result

        # cresult is a sum of Dirac deltas and its derivatives so is known
        # to be causal.

        return cresult, uresult

    def product_undef(self, expr, s, t, **kwargs):

        # Handle expressions with a function of s, e.g., V(s) * Y(s), V(s)
        # / s etc.

        terms = expr.expand().as_ordered_terms()

        result = 0
        for term in terms:
            result += self.product_undef1(term, s, t, **kwargs)
        return result

    def product_undef1(self, expr, s, t, **kwargs):

        const, expr = factor_const(expr, s)

        factors = expr.as_ordered_factors()
        if len(factors) < 2:
            cresult, uresult = self.term1(expr, s, t, **kwargs)
            return const * (cresult + uresult)

        # Help s * 1 / (s + R * C) * I(s)
        if (len(factors) > 2 and not
            isinstance(factors[1], AppliedUndef) and
                isinstance(factors[2], AppliedUndef)):
            factors = [factors[0], factors[2], factors[1]] + factors[3:]

        if isinstance(factors[1], AppliedUndef):
            # Try to expose more simple cases, e.g. (R + s * L) * V(s)
            terms = factors[0].expand().as_ordered_terms()
            if len(terms) >= 2:
                result = Zero
                for term in terms:
                    result += self.product_undef(factors[1]
                                                 * term, s, t, **kwargs)
                return result * const

        cresult, uresult = self.term1(factors[0], s, t, **kwargs)
        result = cresult + uresult

        if kwargs.get('causal', False):
            # Assume that all functions are causal in the expression.
            t1 = Zero
            t2 = t
        else:
            # With unilateral Laplace transform need to set lower
            # limit at 0 rather than -oo.
            t1 = Zero
            t2 = sym.oo

        intnum = 0
        for m in range(len(factors) - 1):
            if m == 0 and isinstance(factors[1], AppliedUndef):
                # Note, as_ordered_factors puts powers of s before the functions.
                if factors[0] == s:
                    # Handle differentiation
                    # Convert s * V(s) to d v(t) / dt
                    result = self.func(factors[1], s, t)
                    result = sym.Derivative(result, t)
                    if not kwargs.get('zero_initial_conditions', True):
                        fname = factors[1].func.__name__
                        func = sym.Function(fname[0].lower() + fname[1:])
                        result += func(0) * sym.DiracDelta(t)
                    continue
                elif factors[0].is_Pow and factors[0].args[0] == s and factors[0].args[1] > 0:
                    # Handle higher order differentiation
                    # Convert s ** 2 * V(s) to d^2 v(t) / dt^2
                    result = self.func(factors[1], s, t)
                    result = sym.Derivative(result, t, factors[0].args[1])
                    if not kwargs.get('zero_initial_conditions', True):
                        fname = factors[1].func.__name__
                        func = sym.Function(fname[0].lower() + fname[1:])
                        v = func(t)
                        order = factors[0].args[1]
                        for m in range(order - 1):
                            result += sym.Derivative(v, t, m).subs(t, 0) * \
                                sym.DiracDelta(t, order - m - 1)
                        result += sym.Derivative(v, t, order - 1).subs(t, 0) * \
                            sym.DiracDelta(t)

                    continue
                elif factors[0].is_Pow and factors[0].args[0] == s \
                        and factors[0].args[1] == -1:
                    # Handle integration  1 / s * V(s)
                    tau = self.dummy_var(expr, 'tau', level=intnum, real=True)
                    intnum += 1
                    result = self.func(factors[1], s, tau)
                    result = sym.Integral(result, (tau, t1, t))
                    continue

            # Convert product to convolution
            tau = self.dummy_var(expr, 'tau', level=intnum, real=True)
            intnum += 1
            cresult, uresult = self.term1(factors[m + 1], s, t, **kwargs)
            expr2 = cresult + uresult
            result = sym.Integral(result.subs(t, t - tau) * expr2.subs(t, tau),
                                  (tau, t1, t2))

        return result * const

    def power(self, expr, s, t):

        # Handle expressions with a power of s.
        if not (expr.is_Pow and expr.args[0] == s):
            self.error('Expression not a power of s')
        exponent = expr.args[1]

        # Have many possible forms; the common ones are:
        # s**a, s**-a, s**(1+a), s**(1-a), s**-(1+a), s**(a-1)
        # Cannot tell if 1-a is positive.

        if exponent.is_positive:
            # Unfortunately, SymPy does not seem to support fractional
            # derivatives...
            return sym.Derivative(sym.DiracDelta(t), t, exponent,
                                  evaluate=False)

        if exponent.is_negative:
            return sym.Pow(t, -exponent - 1) / sym.Gamma(-exponent)

        self.error('Cannot determine sign of exponent')

    def delay_factor(self, expr, var):

        delay = Zero
        rest = One

        for f in expr.as_ordered_factors():
            b, e = f.as_base_exp()
            if b == sym.E and e.is_polynomial(var):
                p = sym.Poly(e, var)
                c = p.all_coeffs()
                if p.degree() == 1:
                    delay -= c[0]
                    if c[1] != 0:
                        rest *= sym.exp(c[1])
                    continue

            rest *= f
        return rest, delay

    def sympy(self, expr, s, t):

        self.debug('Resorting to SymPy')

        # This barfs when needing to generate Dirac deltas
        from sympy.integrals.transforms import inverse_laplace_transform
        result = inverse_laplace_transform(expr, s, t)

        if result.has(sym.InverseLaplaceTransform):
            self.error('SymPy does not know either')
        return result

    def tline_end(self, expr, s, t):
        """Attempt to find response at end of unterminated lossless
        transmission line."""

        # Look for expression of form 1 / (a * cosh(s * T) + b * sinh(s * T))

        func = sym.DiracDelta
        if (expr.is_Mul and expr.args[0].is_Pow and
                expr.args[0].args[0] == s and expr.args[0].args[1] == -1):
            expr *= s
            func = sym.Heaviside

        if not expr.is_Pow or expr.args[1] != -1:
            self.error('Not 1 / X')

        denom = expr.args[0]

        arg = None
        a = Zero
        b = Zero

        for term in denom.expand().as_ordered_terms():

            const, e = factor_const(term, s)
            if e.is_Function and e.func == sym.cosh:
                if arg is None:
                    arg = e.args[0]
                elif arg != e.args[0]:
                    self.error('Mismatch cosh arg')
                a += const
            elif e.is_Function and e.func == sym.sinh:
                if arg is None:
                    arg = e.args[0]
                elif arg != e.args[0]:
                    self.error('Mismatch sinh arg')
                b += const
            else:
                self.error('Term does not include cosh or sinh')

        T = arg / s
        if T.has(s):
            self.error('Frequency dependent delay')

        m = self.dummy_var(expr, 'm', level=0, real=True)
        g = (b - a) / (b + a)
        d = (b + a) / 2
        if g == 0:
            h = 1 / d * func(t - T)
        else:
            h = 1 / d * sym.Sum(g**m * func(t - (2 * m + 1)
                                * T), (m, 0, sym.oo))
        return h

    def tline_start(self, expr, s, t):
        """Attempt to find response at start of unterminated lossless
        transmission line."""

        # Look for expression of form
        # c * cosh(s * T) + d * sinh(s * T) / (a * cosh(s * T) + b * sinh(s * T))
        func = sym.DiracDelta
        if (expr.is_Mul and expr.args[0].is_Pow and
                expr.args[0].args[0] == s and expr.args[0].args[1] == -1):
            expr *= s
            func = sym.Heaviside

        numer, denom = expr.as_numer_denom()

        arg = None
        a = Zero
        b = Zero
        c = Zero
        d = Zero

        for term in denom.expand().as_ordered_terms():

            const, e = factor_const(term, s)
            if e.is_Function and e.func == sym.cosh:
                if arg is None:
                    arg = e.args[0]
                elif arg != e.args[0]:
                    self.error('Mismatch cosh arg')
                a += const
            elif e.is_Function and e.func == sym.sinh:
                if arg is None:
                    arg = e.args[0]
                elif arg != e.args[0]:
                    self.error('Mismatch sinh arg')
                b += const
            else:
                self.error('Term does not include cosh or sinh')

        for term in numer.expand().as_ordered_terms():

            const, e = factor_const(term, s)
            if e.is_Function and e.func == sym.cosh:
                if arg is None:
                    arg = e.args[0]
                elif arg != e.args[0]:
                    self.error('Mismatch cosh arg')
                c += const
            elif e.is_Function and e.func == sym.sinh:
                if arg is None:
                    arg = e.args[0]
                elif arg != e.args[0]:
                    self.error('Mismatch sinh arg')
                d += const
            else:
                self.error('Term does not include cosh or sinh')

        T = arg / s
        if T.has(s):
            self.error('Frequency dependent delay')

        Z0 = sym.sqrt((b * d - a * c) * d / (d**2 - c**2))
        Zl = Z0 * c / d
        Zs = (a - Z0 * Zl) / Z0
        K = d / Z0**2
        # The simplify can be slow
        Gammas = ((Zs - Z0) / (Zs + Z0)).simplify()
        Gammal = ((Zl - Z0) / (Zl + Z0)).simplify()

        m = self.dummy_var(expr, 'm', level=0, real=True)
        h1 = (1 - Gammas) / 2 * func(t)
        if Gammal == 0:
            h2 = 0
        elif Gammas == 0:
            h2 = (1 - Gammas) * (1 + Gammas) * func(t - 2 * m * T)
        else:
            h2 = (1 - Gammas) * (1 + Gammas) / Gammas * \
                sym.Sum((Gammas * Gammal)**m *
                        func(t - 2 * m * T), (m, 1, sym.oo))
        return K * (h1 + h2)

    def term1(self, expr, s, t, **kwargs):

        const, expr = factor_const(expr, s)

        if isinstance(expr, AppliedUndef):
            # Handle V(s), 3 * V(s) etc.  If causal is True it is assumed
            # that the unknown functions are causal.
            result = self.func(expr, s, t)
            return result * const, Zero

        if expr.has(AppliedUndef):
            return const * self.product_undef(expr, s, t, **kwargs), Zero

        try:
            # This is the common case.
            cresult, uresult = self.ratfun(expr, s, t, **kwargs)
            return const * cresult, const * uresult
        except:
            pass

        if expr.is_Pow and expr.args[1] == -1 and expr.args[0].is_Function:
            arg = expr.args[0].args[0]

            m = self.dummy_var(expr, 'm', level=0, real=True)
            if expr.args[0].func == sym.cosh:
                scale, shift = scale_shift(arg, s)
                if shift == 0:
                    return const * 2 * sym.Sum((-1)**m * sym.DiracDelta(t - scale * (2 * m + 1)), (m, 0, sym.oo)), Zero
            elif expr.args[0].func == sym.sinh:
                scale, shift = scale_shift(arg, s)
                if shift == 0:
                    return const * 2 * sym.Sum(sym.DiracDelta(t - scale * (2 * m + 1)), (m, 0, sym.oo)), Zero
            elif expr.args[0].func == sym.tanh:
                scale, shift = scale_shift(arg, s)
                if shift == 0:
                    return const * 2 * sym.Sum((-1)**m * sym.DiracDelta(t - scale * (2 * m + 1)), (m, 1, sym.oo)) + const * sym.DiracDelta(t), Zero

        if expr.has(sym.cosh) and expr.has(sym.sinh):
            try:
                return const * self.tline_end(expr, s, t), Zero
            except:
                pass
            try:
                return const * self.tline_start(expr, s, t), Zero
            except:
                pass

        if expr.is_Pow and expr.args[0] == s:
            return Zero, const * self.power(expr, s, t)

        self.error('Cannot determine inverse Laplace transform')

    def term(self, expr, s, t, **kwargs):

        expr, delay = self.delay_factor(expr, s)

        if delay == 0 and expr.has(sym.exp):
            # Handle cases like 1 / (s**2 * exp(5 * s) + s * exp(5 * s))
            expr1 = expr.simplify()
            expr2, delay2 = self.delay_factor(expr1, s)
            if not expr2.has(sym.exp):
                # Simplify can make things worse, e.g., 1 - exp(-5 *s)
                # becomes exp(-5 * s) * (exp(5 * s) - 1)
                expr = expr2
                delay = delay2

        try:
            cresult, uresult = self.term1(expr, s, t, **kwargs)
        except:

            # With deep=True, SymPy makes mess of (1 - exp(-s * T))
            terms = expr.expand(deep=False).as_ordered_terms()

            if len(terms) > 1:

                try:
                    # See if can convert to convolutions...
                    return self.product_undef(expr, s, t, **kwargs), Zero
                except:
                    pass

                uresult = Zero
                cresult = Zero

                for term in terms:
                    term = term.simplify()
                    cterm, uterm = self.term(term, s, t, **kwargs)
                    cresult += cterm
                    uresult += uterm
                return cresult, uresult

            expr = expr.simplify()

            try:
                cresult, uresult = self.term1(expr, s, t, **kwargs)
            except:
                return Zero, self.sympy(expr, s, t)

        if delay != 0:
            cresult = cresult.subs(t, t - delay)
            uresult = uresult.subs(t, t - delay)

            # h(t) = g(t - T)
            # If h(t) is known to be causal and T >= 0, then g(t)
            # is also causal.  If T > 0, then causality is violated.

            # If h(t) is not known to be causal then we can only infer
            # it for t >= 0 from its Laplace transform H(s).
            # From the delay theorem, H(s) = G(s) * exp(-s * T).
            # So given G(s) we can only infer g(t) for t >= 0.
            # This implies that we can only infer h(t) for t >= T.
            # This creates a can of worms since different terms of
            # an expression can have different conditions, say
            # exp(-3 * s) / s**2 + exp(-4 * s) / s.  So for now,
            # assume causal expression...

            if not kwargs.get('causal', False) and expr.has(AppliedUndef):
                warn('Assuming causal expression')

            if not delay.is_negative:
                if not delay.is_positive:
                    warn('Assuming %s is positive' % delay)
                cresult += uresult * sym.Heaviside(t - delay)
                uresult = Zero
            else:
                self.error('Causality violated with time advance %s.' % delay)

        else:
            if kwargs.get('causal', False):
                cresult += uresult * sym.Heaviside(t)
                uresult = Zero

        return cresult, uresult


inverse_laplace_transformer = InverseLaplaceTransformer()


def inverse_laplace_transform(expr, s, t, evaluate=True, **kwargs):
    """Calculate inverse Laplace of X(s) and return x(t).

    The unilateral Laplace transform cannot determine x(t) for n < 0
    unless given additional information:

    `dc` -- x(t) = constant
    `causal` -- x(t) = 0 for n < 0.
    `ac` -- x(t) = A cos(a * n) + B * sin(b * n)
    """

    return inverse_laplace_transformer.transform(expr, s, t,
                                                 evaluate=evaluate,
                                                 **kwargs)


def ILT(expr, s, t, evaluate=True, **kwargs):
    """Calculate inverse Laplace of X(s) and return x(t).

    The unilateral Laplace transform cannot determine x(t) for n < 0
    unless given additional information:

    `dc` -- x(t) = constant
    `causal` -- x(t) = 0 for n < 0.
    `ac` -- x(t) = A cos(a * n) + B * sin(b * n)
    """

    return inverse_laplace_transformer.transform(expr, s, t,
                                                 evaluate=evaluate,
                                                 **kwargs)

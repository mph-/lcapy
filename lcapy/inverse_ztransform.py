"""This module provides support for the inverse z-transform.  It
calculates the unilateral inverse z-transform.

Copyright 2021--2022 Michael Hayes, UCECE

"""

from .transformer import UnilateralInverseTransformer
from .ratfun import Ratfun
from .utils import factor_const, scale_shift, pair_conjugates
from .sym import sympify, simplify, miscsymbol, AppliedUndef
from .utils import factor_const, scale_shift
from .extrafunctions import UnitImpulse, UnitStep
import sympy as sym
from sympy.simplify.fu import TR6, TR9
from warnings import warn


__all__ = ('IZT', 'inverse_ztransform')


class InverseZTransformer(UnilateralInverseTransformer):

    name = 'inverse z-transform'

    def noevaluate(self, expr, z, n):

        # Construct Cauchy integral.
        self.error('TODO')

    def check(self, expr, z, n, **kwargs):

        if expr.has(n):
            self.error('Expression depends on n')

    def key(self, expr, z, n, **kwargs):
        return (expr, z, n,
                kwargs.get('causal', False),
                kwargs.get('pairs', True),
                kwargs.get('damping', None))

    def func(self, expr, n, z, inverse=False):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], n)

        # Convert v[n] to V(z), etc.
        name = expr.func.__name__
        if inverse:
            func = sym.Function(name[0].lower() + name[1:])
        else:
            func = sym.Function(name[0].upper() + name[1:])

        if not scale.is_constant():
            self.error('Cannot determine if time-expansion or decimation')

        if scale == 1:
            result = func(z)

            if shift != 0:
                result = result * z ** shift
            return result

        if scale.is_integer:
            self.error('Cannot do decimation yet')

        if not scale.is_rational:
            self.error('Cannot handle arbitrary scaling')

        if scale.p != 1:
            self.error('Cannot handle non-integer time-expansion')

        result = func(z ** scale.q)

        if shift != 0:
            result = result * z ** shift
        return result

    def ratfun(self, expr, z, n, **kwargs):

        expr = expr / z

        # Handle special case 1 / (z**m * (z - 1)) since this becomes u[n - m]
        # The default method produces u[n] - delta[n] for u[n-1].  This is correct
        # but can be simplified.
        # In general, 1 / (z**m * (z - a)) becomes a**n * u[n - m]

        if (len(expr.args) == 2 and expr.args[1].is_Pow and
            expr.args[1].args[0].is_Add and
            expr.args[1].args[0].args[0] == -1 and
                expr.args[1].args[0].args[1] == z):

            delay = None
            if expr.args[0] == z:
                delay = 1
            elif expr.args[0].is_Pow and expr.args[0].args[0] == z:
                a = expr.args[0].args[1]
                if a.is_positive:
                    warn('Dodgy z-transform 1.  Have advance of unit step.')
                elif not a.is_negative:
                    warn('Dodgy z-transform 2.  May have advance of unit step.')
                delay = -a
            elif (expr.args[0].is_Pow and expr.args[0].args[0].is_Pow and
                  expr.args[0].args[0].args[0] == z and
                  expr.args[0].args[0].args[1] == -1):
                a = expr.args[0].args[1]
                if a.is_negative:
                    warn('Dodgy z-transform 3.  Have advance of unit step.')
                elif not a.is_positive:
                    warn('Dodgy z-transform 4.  May have advance of unit step.')
                delay = a

            if delay is not None:
                return UnitStep(n - delay), sym.S.Zero

        zexpr = Ratfun(expr, z)

        Q, M, D, delay, undef = zexpr.as_QMA()

        cresult = sym.S.Zero
        uresult = sym.S.Zero

        if Q:
            Qpoly = sym.Poly(Q, z)
            C = Qpoly.all_coeffs()
            for m, c in enumerate(C):
                cresult += c * UnitImpulse(n - len(C) + m + 1)

        # There is problem with determining residues if
        # have 1/(z*(-a/z + 1)) instead of 1/(-a + z).  Hopefully,
        # simplify will fix things...
        expr = (M / D).simplify()
        # M and D may contain common factors before simplification, so redefine M and D
        M = sym.numer(expr)
        D = sym.denom(expr)
        for factor in expr.as_ordered_factors():
            if factor == sym.oo:
                return factor, factor

        zexpr = Ratfun(expr, z, **kwargs)
        poles = zexpr.poles(damping=kwargs.get('damping', None))
        poles_dict = {}
        for pole in poles:
            # Replace cos()**2-1 by sin()**2
            pole.expr = TR6(sym.expand(pole.expr))
            pole.expr = sym.simplify(pole.expr)
            # Remove abs value from sin()
            pole.expr = pole.expr.subs(sym.Abs, sym.Id)
            poles_dict[pole.expr] = pole.n

        # Juergen Weizenecker HsKa

        # Make two dictionaries in order to handle them differently and make
        # pretty expressions
        if kwargs.get('pairs', True):
            pole_pair_dict, pole_single_dict = pair_conjugates(poles_dict)
        else:
            pole_pair_dict, pole_single_dict = {}, poles_dict

        # Make n (=number of poles) different denominators to speed up
        # calculation and avoid sym.limit.  The different denominators are
        # due to shortening of poles after multiplying with (z-z1)**o
        if not (M.is_polynomial(z) and D.is_polynomial(z)):
            print("Numerator or denominator may contain 1/z terms: ", M, D)

        n_poles = len(poles)
        # Leading coefficient of denominator polynom
        a_0 = sym.LC(D, z)
        # The canceled denominator (for each (z-p)**o)
        shorten_denom = {}
        for i in range(n_poles):
            shorten_term = sym.prod(
                [(z - poles[j].expr)**(poles[j].n) for j in range(n_poles) if j != i], a_0)
            shorten_denom[poles[i].expr] = shorten_term

        # Run through single poles real or complex, order 1 or higher
        for pole in pole_single_dict:

            p = pole

            # Number of occurrences of the pole.
            o = pole_single_dict[pole]

            # X(z)/z*(z-p)**o after shortening.
            expr2 = M / shorten_denom[p]

            if o == 0:
                continue

            if o == 1:
                r = sym.simplify(sym.expand(expr2.subs(z, p)))

                if p == 0:
                    cresult += r * UnitImpulse(n)
                else:
                    uresult += r * p ** n
                continue

            # Handle repeated poles.
            all_derivatives = [expr2]
            for i in range(1, o):
                all_derivatives += [sym.diff(all_derivatives[i - 1], z)]

            bino = 1
            sum_p = 0
            for i in range(1, o + 1):
                m = o - i
                derivative = all_derivatives[m]
                # Derivative at z=p
                derivative = sym.expand(derivative.subs(z, p))
                r = sym.simplify(derivative) / sym.factorial(m)

                if p == 0:
                    cresult += r * UnitImpulse(n - i + 1)
                else:
                    sum_p += r * bino * p**(1 - i) / sym.factorial(i - 1)
                    bino *= n - i + 1

            uresult += sym.simplify(sum_p * p**n)

        # Run through complex pole pairs
        for pole in pole_pair_dict:

            p1 = pole[0]
            p2 = pole[1]

            # Number of occurrences of the pole pair
            o1 = pole_pair_dict[pole]
            # X(z)/z*(z-p)**o after shortening
            expr_1 = M / shorten_denom[p1]
            expr_2 = M / shorten_denom[p2]

            # Oscillation parameter
            lam = sym.sqrt(sym.simplify(p1 * p2))
            p1_n = sym.simplify(p1 / lam)
            # term is of form exp(j*arg())
            if len(p1_n.args) == 1 and p1_n.is_Function and p1_n.func == sym.exp:
                omega_0 = sym.im(p1_n.args[0])
            # term is of form cos() + j sin()
            elif p1_n.is_Add and sym.re(p1_n).is_Function and sym.re(p1_n).func == sym.cos:
                p1_n = p1_n.rewrite(sym.exp)
                omega_0 = sym.im(p1_n.args[0])
            # general form
            else:
                omega_0 = sym.simplify(sym.arg(p1_n))

            if o1 == 1:
                r1 = expr_1.subs(z, p1)
                r2 = expr_2.subs(z, p2)

                r1_re = sym.re(r1).simplify()
                r1_im = sym.im(r1).simplify()

                # if pole pairs is selected, r1=r2*

                # Handle real part
                uresult += 2 * TR9(r1_re) * lam ** n * sym.cos(omega_0 * n)
                uresult -= 2 * TR9(r1_im) * lam ** n * sym.sin(omega_0 * n)

            else:
                bino = 1
                sum_b = 0
                # Compute first all derivatives needed
                all_derivatives_1 = [expr_1]
                for i in range(1, o1):
                    all_derivatives_1 += [
                        sym.diff(all_derivatives_1[i - 1], z)]

                # Loop through the binomial series
                for i in range(1, o1 + 1):
                    m = o1 - i

                    # m th derivative at z=p1
                    derivative = all_derivatives_1[m]
                    r1 = derivative.subs(z, p1) / sym.factorial(m)
                    # prefactors
                    prefac = bino * lam ** (1 - i) / sym.factorial(i - 1)
                    # simplify r1
                    r1 = r1.rewrite(sym.exp).simplify()
                    # sum
                    sum_b += prefac * r1 * sym.exp(sym.I * omega_0 * (1 - i))
                    # binomial coefficient
                    bino *= n - i + 1

                # take result = lam**n * (sum_b*sum_b*exp(j*omega_0*n) + cc)
                aa = sym.simplify(sym.re(sum_b))
                bb = sym.simplify(sym.im(sum_b))
                uresult += 2 * (aa * sym.cos(omega_0 * n) -
                                bb * sym.sin(omega_0 * n)) * lam**n

        # cresult is a sum of Dirac deltas and its derivatives so is known
        # to be causal.

        return cresult, uresult

    def product(self, expr, z, n, **kwargs):

        # Handle expressions with a function of z, e.g., V(z) * Y(z), V(z)
        # / z etc.

        zsym = sympify(str(z))

        if kwargs.get('causal', False):
            # Assume that all functions are causal in the expression.
            n1 = sym.S.Zero
            n2 = n
        else:
            n1 = -sym.oo
            n2 = sym.oo

        const, expr = factor_const(expr, z)

        factors = expr.as_ordered_factors()
        if len(factors) < 2:
            self.error('Expression does not have multiple factors')

        if (len(factors) == 3 and factors[0] == z and
            isinstance(factors[2], AppliedUndef) and
            factors[1].is_Pow and factors[1].args[1] == -1 and
            factors[1].args[0].is_Add and factors[1].args[0].args[0] == -1
                and factors[1].args[0].args[1] == z):
            # Handle cumulative sum  z / (z - 1) * V(z)
            m = self.dummy_var(expr, 'm', level=0, real=True)
            result = self.func(factors[2], z, m, True)
            result = sym.Sum(result, (m, n1, n))
            return result

        # TODO, is this useful?
        if (len(factors) > 2 and not
            isinstance(factors[1], AppliedUndef) and
                isinstance(factors[2], AppliedUndef)):
            factors = [factors[0], factors[2], factors[1]] + factors[3:]

        # TODO, is this useful?
        if isinstance(factors[1], AppliedUndef):
            terms = factors[0].as_ordered_terms()
            if len(terms) >= 2:
                result = sym.S.Zero
                for term in terms:
                    result += self.product(factors[1] * term, z, n, **kwargs)
                return const * result

        cresult, uresult = self.term1(factors[0], z, n, **kwargs)
        result = cresult + uresult

        intnum = 0
        for m in range(len(factors) - 1):
            if m == 0 and isinstance(factors[1], AppliedUndef):
                # Note, as_ordered_factors puts powers of z before the functions.
                if factors[0] == z:
                    # Handle time-advance
                    # Convert z * V(z) to v[n + 1]
                    # TODO, fix for unilateral ZT
                    result = self.func(factors[1], z, n, True)
                    result = result.subs(n, n + 1)
                    continue
                elif factors[0].is_Pow and factors[0].args[0] == z and factors[0].args[1] > 0:
                    # Handle higher order advances
                    # Convert z ** k * V(z) to v[n + k]
                    result = self.func(factors[1], z, n, True)
                    result = result.subs(n, n + factors[0].args[1])
                    continue
                elif factors[0].is_Pow and factors[0].args[0] == z and factors[0].args[1] < 0:
                    # Handle time-delay  1 / z ** k * V(z)
                    result = self.func(factors[1], z, n, True)
                    result = result.subs(n, n + factors[0].args[1])
                    continue
                elif (factors[0].is_Pow and
                      factors[0].args[0].is_Add and
                      factors[0].args[1] == -1 and
                      factors[0].args[0].args[0] == 1 and
                      factors[0].args[0].args[1].is_Mul and
                      factors[0].args[0].args[1].args[0] == -1 and
                      factors[0].args[0].args[1].args[1].is_Pow and
                      factors[0].args[0].args[1].args[1].args[1] == -1 and
                      factors[0].args[0].args[1].args[1].args[0] is zsym):
                    # Handle cumulative sum  1 / (1 - 1 / z) * V(z)
                    m = self.dummy_var(expr, 'm', level=intnum, real=True)
                    result = self.func(factors[1], z, m, True)
                    intnum += 1
                    result = sym.Sum(result, (m, n1, n))
                    continue

            # Convert product to convolution
            dummy = self.dummy_var(expr, 'm', level=intnum, real=True)
            intnum += 1
            cresult, uresult = self.term1(factors[m + 1], z, n, **kwargs)
            expr2 = cresult + uresult
            kernel = result.subs(n, n - dummy) * expr2.subs(n, dummy)
            sresult = sym.Sum(kernel, (dummy, n1, n2))
            result = sresult

        return const * result

    def power(self, expr, z, n):

        # Handle expressions with a power of z.
        if (expr.is_Pow and expr.args[0] == z):
            exponent = expr.args[1]

            if exponent.is_positive:
                warn('Dodgy z-transform.  Have advance of unit impulse.')
            elif not exponent.is_negative:
                warn('Dodgy z-transform.  May have advance of unit impulse.')

            return UnitImpulse(n + exponent), sym.S.Zero

        # Handle expressions with a power of (1 / z).
        if (expr.is_Pow and expr.args[0].is_Pow and
                expr.args[0].args[0] == z and expr.args[0].args[1] == -1):

            exponent = expr.args[1]

            if exponent.is_negative:
                warn('Dodgy z-transform.  Have advance of unit impulse.')
            elif not exponent.is_positive:
                warn('Dodgy z-transform.  May have advance of unit impulse.')

            return UnitImpulse(n - exponent), sym.S.Zero

        self.error('Expression is not a power of z')

    def term1(self, expr, z, n, **kwargs):

        const, expr = factor_const(expr, z)

        if expr == 1:
            return const * UnitImpulse(n), sym.S.Zero

        if isinstance(expr, AppliedUndef):
            # Handle V(z), 3 * V(z) etc.  If causal is True it is assumed
            # that the unknown functions are causal.
            result = self.func(expr, z, n, True)
            return const * result, sym.S.Zero

        if expr.has(AppliedUndef):
            return const * self.product(expr, z, n, **kwargs), sym.S.Zero

        if expr == z:
            warn('Dodgy z-transform.  Have advance of unit impulse.')
            return const * UnitImpulse(n + 1), sym.S.Zero

        if (expr.is_Pow and
            (expr.args[0] == z or
             (expr.args[0].is_Pow and
              expr.args[0].args[0] == z and expr.args[0].args[1] == -1))):
            cresult, uresult = self.power(expr, z, n)
            return const * cresult, const * uresult

        try:
            # This is the common case.
            cresult, uresult = self.ratfun(expr, z, n)
            return const * cresult, const * uresult
        except:
            pass

        # As last resort see if can convert to convolutions...
        return sym.S.Zero, const * self.product(expr, z, n, **kwargs)

    def term(self, expr, z, n, **kwargs):

        cresult, uresult = self.term1(expr, z, n, **kwargs)

        if kwargs.get('causal', False):
            uresult = uresult * UnitStep(n)

        return cresult, uresult


inverse_ztransformer = InverseZTransformer()


def inverse_ztransform(expr, z, n, evaluate=True, **kwargs):
    """Calculate inverse z-Transform of X(s) and return x[n].

    The unilateral z-Transform transform cannot determine x[n] for n < 0
    unless given additional information:

    `dc` -- x[n] = constant
    `causal` -- x[n] = 0 for n < 0.
    `ac` -- x[n] = A cos(a * n) + B * sin(b * n)
    """

    return inverse_ztransformer.transform(expr, z, n, evaluate=evaluate,
                                          **kwargs)


def IZT(expr, z, n, evaluate=True, **kwargs):
    """Calculate inverse z-Transform of X(s) and return x[n].

    The unilateral z-Transform transform cannot determine x[n] for n < 0
    unless given additional information:

    `dc` -- x[n] = constant
    `causal` -- x[n] = 0 for n < 0.
    `ac` -- x[n] = A cos(a * n) + B * sin(b * n)
    """

    return inverse_ztransformer.transform(expr, z, n, evaluate=evaluate,
                                          **kwargs)

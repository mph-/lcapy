"""
This module provides support for rational functions.

Copyright 2016--2022 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from .sym import sympify, AppliedUndef
from .cache import lru_cache
from .utils import pair_conjugates, factor_const

Zero = sym.S.Zero
One = sym.S.One

class Pole(object):

    def __init__(self, expr, n, damping=None):
        self.damping = damping
        self.n = n
        self.orig_expr = expr
        self.d = self._decompose(expr)
        self.expr = self._rewrite()

    def _decompose_bar(self, expr):

        # Look for scale * sqrt(dexpr) where dexpr = aexpr - bexpr

        dexpr = None
        scale = One

        for factor in expr.as_ordered_factors():
            if (factor.is_Pow and factor.args[1] == One / 2 and
                factor.args[0].is_Add and factor.args[0].args[1].is_Mul and
                factor.args[0].args[1].args[0] < 0):
                if dexpr is not None:
                    return None
                dexpr = factor.args[0]
            else:
                scale *= factor
        if dexpr is None:
            return None
        return scale, dexpr

    def _decompose_foo(self, expr):

        # Look for offset + scale * sqrt(dexpr) where dexpr = aexpr - bexpr

        offset = Zero
        scale = One
        found = False
        dexpr = None

        for term in expr.as_ordered_terms():
            p = self._decompose_bar(term)
            if p is None:
                offset += term
            elif dexpr is not None:
                return None
            else:
                scale, dexpr = p

        return offset, scale, dexpr

    def _decompose(self, expr):

        scale = One
        scale2 = One
        offset = One
        dexpr = None

        for factor in expr.as_ordered_factors():
            if factor.is_Add:
                p = self._decompose_foo(factor)
                if p is None:
                    scale *= factor
                elif dexpr is not None:
                    return None
                else:
                    offset, scale2, dexpr = p
            else:
                scale *= factor

        if dexpr is None:
            dexpr = Zero
        return scale, offset, scale2, dexpr

    @property
    def conj(self):
        return self.conjugate()

    def conjugate(self):
        return self._rewrite(conjugate=True)

    def _rewrite(self, conjugate=False):

        if self.damping in (None, 'over') or self.d is None:
            return self.orig_expr.conjugate() if conjugate else self.orig_expr

        d = self.d

        if self.damping == 'under':
            if conjugate:
                return d[0] * (d[1] - 1j * d[2] * sym.sqrt(-d[3]))
            else:
                return d[0] * (d[1] + 1j * d[2] * sym.sqrt(-d[3]))

        elif self.damping == 'critical':
            # This puts constraints on variables since d[2] == 0.
            return (d[0] * d[1]).conjugate() if conjugate else d[0] * d[1]
        else:
            raise ValueError('Unknown damping %s' % self.damping)


def as_numer_denom_poly(expr, var):

    N = One
    D = One
    for f in expr.as_ordered_factors():
        if f.is_Pow and f.args[1] == -1:
            D *= f.args[0]
        else:
            N *= f

    try:
        Dpoly = sym.Poly(D, var)
        Npoly = sym.Poly(N, var)
    except:
        N, D = expr.as_numer_denom()
        Dpoly = sym.Poly(D, var)
        Npoly = sym.Poly(N, var)

    return Npoly, Dpoly


def as_numer_denom(expr, var):

    if expr.has(1 / var):
        expr = expr.cancel()

    N = One
    D = One
    for f in expr.as_ordered_factors():
        if f.is_Pow and f.args[1] == -1:
            D *= f.args[0]
        else:
            N *= f

    try:
        Dpoly = sym.Poly(D, var)
        Npoly = sym.Poly(N, var)
        return N, D
    except:
        return expr.as_numer_denom()


def _zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists (or dictionaries) of zeros
    and poles, and from a constant gain

    """

    K = sympify(K)
    zeros = sympify(zeros)
    poles = sympify(poles)

    if isinstance(zeros, (tuple, list)):
        zz = [(var - z) for z in zeros]
    else:
        zz = [(var - z) ** zeros[z] for z in zeros]

    if isinstance(zeros, (tuple, list)):
        pp = [1 / (var - p) for p in poles]
    else:
        pp = [1 / (var - p) ** poles[p] for p in poles]

    return sym.Mul(K, *(zz + pp), evaluate=False)


def _tc2tf(zeros, poles, K=1, var=None):
    """Create a transfer function in time-constant form from lists (or
    dictionaries) of zeros and poles, and from a constant gain

    """

    K = sympify(K)
    zeros = sympify(zeros)
    poles = sympify(poles)

    zz = []
    pp = []

    if isinstance(zeros, (tuple, list)):
        for z in zeros:
            if z == 0:
                zz.append(var)
            else:
                zz.append((var / -z + 1))
                K *= z
    else:
        for z, o in zeros.items():
            if z == 0:
                zz.append(var ** o)
            else:
                zz.append((var / -z + 1) ** zeros[z])
                K *= z ** o

    if isinstance(zeros, (tuple, list)):
        for p in poles:
            if p == 0:
                pp.append(p)
            else:
                pp.append(1 / (var / -p + 1))
                K /= p
    else:
        for p, o in poles.items():
            if p == 0:
                pp.append(var ** o)
            else:
                pp.append(1 / (var / -p + 1) ** o)
                K /= p ** o

    K = K.simplify()
    return sym.Mul(K, *(zz + pp), evaluate=False)


def _pr2tf(poles, residues, var=None):
    """Create a transfer function from lists of poles and residues.

    """
    poles = sympify(poles)
    residues = sympify(residues)

    return sym.Add(*[r / (var - p) for r, p in zip(residues, poles)], evaluate=False)


def as_ratfun_delay(expr, var):
    delay = Zero

    if expr.is_rational_function(var):
        N, D = expr.as_numer_denom()
        return N, D, delay

    F = sym.factor(expr).as_ordered_factors()

    rf = One
    for f in F:
        b, e = f.as_base_exp()
        if b == sym.E and e.is_polynomial(var):
            p = sym.Poly(e, var)
            c = p.all_coeffs()
            if p.degree() == 1:
                delay -= c[0]
                if c[1] != 0:
                    rf *= sym.exp(c[1])
                continue

        rf *= f

    if not rf.is_rational_function(var):
        raise ValueError('Expression not a product of rational function'
                         ' and exponential')

    N, D = rf.as_numer_denom()
    return N, D, delay


def as_ratfun_delay_undef(expr, var):
    delay = Zero
    undef = One

    if expr.is_rational_function(var):
        N, D = as_numer_denom(expr, var)
        return N, D, delay, undef

    F = sym.factor(expr).as_ordered_factors()

    rf = One
    for f in F:
        b, e = f.as_base_exp()
        if b == sym.E and e.is_polynomial(var):
            p = sym.Poly(e, var)
            c = p.all_coeffs()
            if p.degree() == 1:
                delay -= c[0]
                if c[1] != 0:
                    rf *= sym.exp(c[1])
                continue
        if isinstance(f, AppliedUndef):
            undef *= f
            continue

        rf *= f

    if not rf.is_rational_function(var):
        raise ValueError('Expression not a product of rational function,'
                         ' exponential, and undefined functions')

    N, D = rf.as_numer_denom()
    return N, D, delay, undef


class Ratfun(object):

    def __init__(self, expr, var):
        # Don't use cancel, it can cause a mess when have exp(s * T).
        self.expr = expr
        self.var = var

    def as_ratfun_delay(self):
        """Split expr as (N, D, delay)
        where expr = (N / D) * exp(var * delay)

        Note, delay only represents a delay when var is s."""

        return as_ratfun_delay(self.expr, self.var)

    def as_ratfun_delay_undef(self):
        """Split expr as (N, D, delay, undef)
        where expr = (N / D) * exp(var * delay) * undef
        and where N is a polynomial in var,
        D is a polynomial in var, and undef is the product
        of undefined functions, e.g., V(s).

        Note, delay only represents a delay when var is s."""

        return as_ratfun_delay_undef(self.expr, self.var)

    def as_const_undef_rest(self):
        """Split expr as (const, undef, rest)
        where expr = const * undef * rest"""

        expr = self.expr
        var = self.var

        const = One
        undef = One
        rest = One

        F = sym.factor(expr).as_ordered_factors()

        for f in F:
            if isinstance(f, AppliedUndef):
                undef *= f
                continue
            if not f.has(var):
                const *= f
                continue
            rest *= f

        return const, undef, rest

    @lru_cache()
    def roots(self):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        return sym.roots(sym.Poly(self.expr, self.var))

    @lru_cache()
    def zeros(self):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        return Ratfun(self.numerator, self.var).roots()

    @lru_cache()
    def poles(self, damping=None):
        """Return poles of expression as a dictionary of Pole objects.
        Note this may not find all the poles."""

        poles = []
        for p, n in Ratfun(self.denominator, self.var).roots().items():

            pole = Pole(p, n=n, damping=damping)
            for q in poles:
                if q.expr == pole.expr:
                    q.n += pole.n
                    pole.n = 0
                    break
            if pole.n != 0:
                poles.append(pole)

        return poles

    def residue(self, pole, poles):
        """Determine residue for given pole."""

        expr = self.expr
        var = self.var

        # Remove occurrence of pole; sym.cancel
        # doesn't always work, for example, for complex poles.
        occurrences = []
        for p in poles:
            occurrences += [p.n - 1 if p.expr == pole else p.n]

        numer, denom = expr.as_numer_denom()
        Dpoly = sym.Poly(denom, var)
        K = Dpoly.LC()

        D = [(var - p.expr) ** o for p, o in zip(poles, occurrences)]
        denom = sym.Mul(K, *D)

        # Could calculate all residues simultaneously using
        # system of linear equations.

        def method1(numer, denom, var, pole):

            d = sym.limit(denom, var, pole)

            if d != 0:
                tmp = (numer / denom).simplify()

                # Sometimes this takes ages...
                return sym.limit(tmp, var, pole)

            # Use l'Hopital's rule
            tmp = numer / denom
            tmp = sym.diff(tmp, var)

            return sym.limit(tmp, var, pole)

        def method2(numer, denom, var, pole):

            d = denom.subs(var, pole)
            n = numer.subs(var, pole)

            if d != 0:
                return n / d

            ddenom = sym.diff(denom, var)
            return n / ddenom.subs(pole)

        m2 = method2(numer, denom, var, pole)
        return m2.cancel()

    @property
    def numerator_denominator(self):
        """Return numerator and denominator of rational function"""

        return as_numer_denom(self.expr, self.var)

    @property
    def N(self):
        return self.numerator

    @property
    def D(self):
        return self.denominator

    @property
    def numerator(self):
        """Return numerator of rational function"""

        numer, denom = self.numerator_denominator
        return numer

    @property
    def denominator(self):
        """Return denominator of rational function"""

        numer, denom = self.numerator_denominator
        return denom

    def canonical(self, factor_const=True):
        """Convert rational function to canonical form; this is like general
        form but with a unity highest power of denominator.  For
        example,

        (5 * s**2 + 5 * s + 5) / (s**2 + 4)

        If factor_const is True, factor constants from numerator, for example,

        5 * (s**2 + s + 1) / (s**2 + 4)

        See also general, partfrac, standard, timeconst, and ZPK

        """

        try:
            N, D, delay, undef = self.as_ratfun_delay_undef()
        except ValueError:
            # TODO: copy?
            return self.expr

        var = self.var
        Dpoly = sym.Poly(D, var)
        Npoly = sym.Poly(N, var)

        if factor_const:
            K = sym.cancel(Npoly.LC() / Dpoly.LC())
            if delay != 0:
                K *= sym.exp(self.var * delay)

            # Divide by leading coefficient
            N = Npoly.monic().as_expr()
            D = Dpoly.monic().as_expr()

            if D == 1:
                expr = N
            else:
                expr = sym.Mul(N, 1 / D, evaluate=False)

            if K != 1:
                expr = sym.Mul(K, expr, evaluate=False)
            expr *= undef
        else:
            C = Dpoly.LC()
            D = Dpoly.monic().as_expr()
            N = (Npoly.as_expr() / C).simplify()
            if D == 1:
                expr = N
            else:
                if N == 1:
                    expr = 1 / D
                else:
                    expr = sym.Mul(N, 1 / D, evaluate=False)
            if delay != 0:
                expr *= sym.exp(self.var * delay)
            expr *= undef

        return expr

    def general(self):
        """Convert rational function to general form.

        See also canonical, partfrac, standard, timeconst, and ZPK"""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        expr = sym.cancel(N / D, self.var)
        if delay != 0:
            expr *= sym.exp(self.var * delay)

        return expr * undef

    def expandcanonical(self):
        """Expand in terms for different powers with each term
        expressed in canonical form."""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        Npoly = sym.Poly(N, self.var)

        expr = Zero

        for m, c in enumerate(reversed(Npoly.all_coeffs())):
            term = sym.Mul(c.simplify() * self.var ** m, 1 / D)
            expr += term

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return expr * undef

    def partfrac(self, combine_conjugates=False, damping=None, split=True):
        """Convert rational function into partial fraction form.

        If combine_conjugates is True then the pair of partial
        fractions for complex conjugate poles are combined.

        See also canonical, standard, general, timeconst, and ZPK

        """

        try:
            Q, R, D, delay, undef = self.as_QRD(combine_conjugates, damping)
        except ValueError:
            if not split:
                raise

            # Try splitting into terms
            result = 0
            for term in self.expr.as_ordered_terms():
                result += Ratfun(term, self.var).partfrac(combine_conjugates,
                                                          split=False)
            return result

        result = Q
        for R, D in zip(R, D):
            result += R / D

        if delay != 0:
            result *= sym.exp(-self.var * delay)

        result *= undef
        return result

    def standard(self, split=True):
        """Convert rational function into mixed fraction form.

        This is the sum of strictly proper rational function and a
        polynomial.

        See also canonical, general, partfrac, timeconst, and ZPK"""

        try:
            Q, M, D, delay, undef = self.as_QMD()
        except ValueError:
            if not split:
                raise

            # Try splitting into terms
            result = 0
            for term in self.expr.as_ordered_terms():
                result += Ratfun(term, self.var).standard(False)
            return result

        expr = Q + sym.cancel(M / D, self.var)

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return expr * undef

    def timeconst(self):
        """Convert rational function to time constant form with unity
        lowest power of denominator.

        See also canonical, general, partfrac, standard, and ZPK"""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        var = self.var
        Npoly = sym.Poly(N, var)
        Dpoly = sym.Poly(D, var)

        K = Dpoly.EC()

        D = D / K
        N = N / K
        return sym.Mul(N, sym.Pow(D, -1), evaluate=False) * sym.exp(self.var * delay) * undef

    def ZPK(self, combine_conjugates=False):
        """Convert to zero-pole-gain (ZPK) form.

        See also canonical, general, standard, timeconst, and partfrac"""

        zeros, poles, K, undef = self.as_ZPK()

        var = self.var

        if not combine_conjugates:
            return _zp2tf(zeros, poles, K, var) * undef

        pole_pairs, pole_singles = pair_conjugates(poles)
        zero_pairs, zero_singles = pair_conjugates(zeros)

        result1 = 1
        num = 1
        for zeros, order in zero_pairs.items():
            for m in range(order):
                num *= (var**2 - zeros[0] * var - zeros[1] * var + zeros[0] * zeros[1]).simplify()

        den = 1
        for poles, order in pole_pairs.items():
            for m in range(order):
                den *= (var**2 - poles[0] * var - poles[1] * var + poles[0] * poles[1]).simplify()

        result1 *= (num / den)

        result2 = _zp2tf(zero_singles, pole_singles, 1, var) * undef
        result = K * result1 * result2

        return result

    def residues(self, combine_conjugates=False, damping=None):
        """Return residues of partial fraction expansion.

        This is not much use without the corresponding poles.
        It is better to use as_QRD."""

        Q, R, D, delay, undef = self.as_QRD(combine_conjugates, damping)
        return R

    def coeffs(self):

        N, D, delay, undef = self.as_ratfun_delay_undef()

        var = self.var
        Npoly = sym.Poly(N, var)
        Dpoly = sym.Poly(D, var)

        return Npoly.all_coeffs(), Dpoly.all_coeffs()

    @property
    def degree(self):
        """Return the degree (order) of the rational function.

        This the maximum of the numerator and denominator degrees.
        Note zero has a degree of -inf."""

        Npoly, Dpoly = as_numer_denom_poly(self.expr, self.var)
        return max(Npoly.degree(), Dpoly.degree())

    @property
    def Ndegree(self):
        """Return the degree (order) of the numerator of a rational function.
        Note zero has a degree of -inf."""

        Npoly, Dpoly = as_numer_denom_poly(self.expr, self.var)
        return Npoly.degree()

    @property
    def Ddegree(self):
        """Return the degree (order) of the denominator of a rational function.
        Note zero has a degree of -inf."""

        Npoly, Dpoly = as_numer_denom_poly(self.expr, self.var)
        return Dpoly.degree()

    @property
    def is_strictly_proper(self):
        """Return True if the degree of the dominator is greater
        than the degree of the numerator."""

        return self.Ddegree > self.Ndegree

    def as_ZPK(self):
        """Decompose expression into zeros, poles, gain, undef where

        expression = K * (prod_n (var - z_n) / (prod_n (var - p_n)) * undef
        """

        N, D, delay, undef = self.as_ratfun_delay_undef()

        var = self.var
        Npoly = sym.Poly(N, var)
        Dpoly = sym.Poly(D, var)

        K = sym.cancel(Npoly.LC() / Dpoly.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        zeros = sym.roots(Npoly)
        poles = sym.roots(Dpoly)

        return zeros, poles, K, undef

    def as_QMD(self):
        """Decompose expression into Q, M, D, delay, undef where

        `expression = (Q + M / D) * exp(-delay * var) * undef`"""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        # Perform polynomial long division so expr = Q + M / D
        Q, M = sym.div(N, D, self.var)

        return Q, M, D, delay, undef

    def as_QRD_old(self, combine_conjugates=False, damping=None):
        """Decompose expression into Q, R, D, delay, undef where

        expression = (Q + sum_n R_n / D_n) * exp(-delay * var) * undef"""

        Q, M, D, delay, undef = self.as_QMD()

        expr = M / D
        var = self.var

        sexpr = Ratfun(expr, var)
        poles = sexpr.poles(damping=damping)

        if damping == 'critical':
            # Critical damping puts constraints on variables.  If
            # these variables do not occur in numerator of expr
            # perhaps the following code will work.
            raise ValueError('Critical damping not supported')

        polesdict = {}
        for pole in poles:
            polesdict[pole.expr] = pole.n

        R = []
        D = []

        for pole in poles:

            # Number of occurrences of the pole.
            o = polesdict[pole.expr]
            if o == 0:
                continue

            p = pole.expr
            pc = pole.conjugate
            if combine_conjugates and pc != p and pc in polesdict:
                polesdict[pc] -= 1

                D2 = sym.simplify(var**2 - (p + pc) * var + p * pc)

                if o == 1:
                    r = sexpr.residue(p, poles)
                    rc = sexpr.residue(pc, poles)

                    r = sym.simplify(r * (var - pc) + rc * (var - p))
                    R.append(r)
                    D.append(D2)
                else:
                    # Handle repeated complex pole pairs.
                    expr2 = expr * (var - p) ** o
                    for n in range(1, o + 1):
                        m = o - n

                        dexpr = sym.diff(expr2, var, m)
                        r = sym.limit(dexpr, var, p) / sym.factorial(m)

                        rc = r.conjugate()
                        r = sym.simplify(r * (var - pc) ** n + rc * (var - p) ** n)
                        R.append(r)
                        D.append(D2 ** n)
            else:
                # SymPy reorders the following, for example, s - -1 => 1 + s
                # D2 = var - p
                if p == 0:
                    D2 = var
                else:
                    D2 = sym.Add(var, -p, evaluate=False)

                if o == 1:
                    r = sexpr.residue(p, poles)
                    R.append(r)
                    D.append(D2)
                else:
                    # Handle repeated real poles.
                    expr2 = expr * (var - p) ** o
                    for n in range(1, o + 1):
                        m = o - n

                        dexpr = sym.diff(expr2, var, m).simplify()
                        r = sym.limit(dexpr, var, p) / sym.factorial(m)

                        R.append(r)
                        D.append(D2 ** n)

        return Q, R, D, delay, undef

    def _prune_zero_residues(self, R, P, O):

        Rnew = []
        Pnew = []
        Onew = []
        for r, p, o in zip(R, P, O):
            if r != 0:
                Rnew.append(r)
                Pnew.append(p)
                Onew.append(o)
        return Rnew, Pnew, Onew

    def as_QRPO(self, damping=None):
        """Decompose expression into Q, R, P, O, delay, undef where

        `expression = (Q + sum_n R_n / (var - P_n)**O_n) * exp(-delay * var) * undef`
        """

        from .matrix import matrix_inverse

        Q, M, D, delay, undef = self.as_QMD()

        expr = M / D
        var = self.var

        sexpr = Ratfun(expr, var)
        poles = sexpr.poles(damping=damping)

        if damping == 'critical':
            # Critical damping puts constraints on variables.  If
            # these variables do not occur in numerator of expr
            # perhaps the following code will work.
            raise ValueError('Critical damping not supported')

        if len(poles) == 0:
            return Q, [], [], [], delay, undef

        elif len(poles) == 1 and poles[0].n == 1:
            p = poles[0].expr
            d = var - p
            r = (expr * d).cancel()
            return Q, [r], [p], [1], delay, undef

        # Find residues by solving system of equations (equating coefficients
        # method).  For example,  consider an expression with repeated poles:
        # E(s) = N(S) / ((s - p1)**2 * (s - p2))
        #
        # This can be expressed as partial fractions:
        # E(s) = R_1 / (s - p1)**2 + R_2 / (s - p1) + R_3 / (s - p2)
        # where R_n are the residues we would like to find.
        #
        # The factored denominator expression is
        # D(s) = ((s - p1)**2 * (s - p2))
        # and so the numerator expression is
        # N(s) = E(s) * D(s)
        #
        # Expanding the partial fractions to have a common denominator gives:
        # N(s) = R_1 * (s - p2) + R_2 * (s - p1) * (s - p2) + R * (s - p1)**2
        #
        # Equating N(s) with E(s) * D(s) and equating powers of s gives
        # system of equation to find the residues.

        syms = []
        D = []
        O = []
        P = []
        i = 1
        denom_factored = One
        for pole in poles:
            for m in range(pole.n):
                A = sym.Symbol('A_%d' % i)
                p = pole.expr
                d = (var - p) ** (m + 1)
                syms.append(A)
                D.append(d)
                P.append(p)
                O.append(m + 1)
                i += 1
            denom_factored *= (var - pole.expr) ** pole.n

        rhs = Zero
        for denom, residue in zip(D, syms):
            # Could be more cunning to avoid using cancel
            rhs += residue * (denom_factored / denom).cancel()

        # Note, the following can differ from self.numerator by
        # a scale factor.  TODO, find a better way to determine
        # the scale factor.
        lhs = (expr * denom_factored).cancel()

        rhspoly = sym.poly(rhs, var)
        lhspoly = sym.poly(lhs, var)

        rc = rhspoly.all_coeffs()
        lc = lhspoly.all_coeffs()
        if len(lc) < len(rc):
            lc = [0] * (len(rc) - len(lc)) + lc

        A, _ = sym.linear_eq_to_matrix(rc, syms)

        # Solve system of equations to find residues.
        R = list(matrix_inverse(A) * sym.Matrix(lc))

        # Remove elements where the residue is zero.
        R, P, O = self._prune_zero_residues(R, P, O)

        return Q, R, P, O, delay, undef

    def _combine_conjugates(self, R, D, P, O):

        Rnew = []
        Dnew = []

        for m in range(len(R)):

            r = R[m]
            if r == 0:
                continue
            p = P[m]
            o = O[m]
            rc = r.conjugate()
            pc = p.conjugate()

            has_conjugate = False
            for n in range(m + 1, len(R)):
                if o != O[n] or pc != P[n] or rc != R[n]:
                    continue
                if o > 1:
                    continue
                R[n] = 0
                has_conjugate = True
                break

            if has_conjugate:
                Rnew.append(((self.var - pc) * r + (self.var - p) * rc).expand())
                Dnew.append(((self.var - p) * (self.var - pc)).expand())
            else:
                Rnew.append(r)
                Dnew.append(D[m])
        return Rnew, Dnew


    def as_QRD(self, combine_conjugates=False, damping=None):
        """Decompose expression into Q, R, D, delay, undef where

        expression = (Q + sum_n R_n / D_n) * exp(-delay * var) * undef"""

        Q, R, P, O, delay, undef = self.as_QRPO(damping)

        D = [(self.var - p)**o for p, o in zip(P, O)]

        if combine_conjugates:
            R, D = self._combine_conjugates(R, D, P, O)

        return Q, R, D, delay, undef

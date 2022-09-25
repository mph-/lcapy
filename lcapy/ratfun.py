"""
This module provides support for rational functions.

Copyright 2016--2022 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from .sym import sympify, AppliedUndef
from .cache import lru_cache, cached_property
from .utils import pair_conjugates, factor_const
from warnings import warn

Zero = sym.S.Zero
One = sym.S.One

# Perhaps have 3 classes to represent responses?
# 1. RationalFunction            N(s) / D(s)
# 2. RationalFunctionWithExp     B(s) * exp(-s * T) / D(s)
# 3. RationalFunctionWithExpSum  Sum(B_m(s) * exp(-s * T_m), m) / D(s)
#
# 2 is a generalized form of 1 with T = 0.


def polyroots(poly, var):
    """Return roots of polynomial `poly` for variable `var`."""

    roots = sym.roots(poly)
    num_roots = 0
    for root, n in roots.items():
        num_roots += n
    if num_roots != poly.degree():
        if poly.degree() >= 5:
            warn('Cannot find symbolic roots for polynomial of degree %d, see Abel-Ruffini theorem' % poly.degree())

        # If the coefficients of the polynomial are numerical,
        # the SymPy nroots function can be used to find
        # numerical approximations to the roots.
        a = set()
        a.add(var)
        if poly.free_symbols == a:
            warn('Only %d of %d roots found, using numerical approximation' %
                 (num_roots, poly.degree()))

            nroots = poly.nroots()

            roots = {}
            for root in nroots:
                if root in roots:
                    roots[root] += 1
                else:
                    roots[root] = 1

            return roots
        warn('Only %d of %d roots found' % (num_roots, poly.degree()))

    return roots


def roots(expr, var):
    """Return roots of expression `expr` for variable `var`."""

    return polyroots(sym.Poly(expr, var), var)


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


def as_B_A_delay_undef(expr, var):

    delay = Zero
    undef = One

    # This does not detect expressions of the form:
    # A(s) * (B(s) - exp(-s * T) / D(s)

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
        """This represents a generalized rational function of the form:
        `N(var) / D(var) = (B(var) / A(var)) * exp(-var * delay) * U(var)`

        where B and D are polynomials in var and U is the product of
        undefined functions of var.

        Note, delay only represents a delay when var is s.

        It does not detect expressions such as:
        `B(var) * (1 - exp(-var * delay)) / A(var)

        """

        self.expr = expr
        self.var = var

        self.B, self.A, self.delay, self.undef = as_B_A_delay_undef(expr, var)

        N = self.B * self.undef
        if self.delay != 0:
            N *= sym.exp(-self.delay * var)

        # These may change...  Perhaps, A should not have a constant factor?
        self.N = N
        self.D = self.A

    def as_B_A_delay_undef(self):

        return self.B, self.A, self.delay, self.undef

    def _roots(self, poly):

        return polyroots(poly, self.var)

    @lru_cache()
    def roots(self):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        return self._roots(sym.Poly(self.expr, self.var))

    @lru_cache()
    def zeros(self):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        return self._roots(self.Bpoly)

    @lru_cache()
    def poles(self, damping=None):
        """Return poles of expression as a list of Pole objects.
        Note this may not find all the poles."""

        roots = self._roots(self.Apoly)

        poles = []

        for p, n in roots.items():

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

    @cached_property
    def Bpoly(self):
        return sym.Poly(self.B, self.var)

    @cached_property
    def Apoly(self):
        return sym.Poly(self.A, self.var)

    def canonical(self, factor_const=True):
        """Convert rational function to canonical form; this is like general
        form but with a unity highest power of denominator.  For
        example,

        (5 * s**2 + 5 * s + 5) / (s**2 + 4)

        If factor_const is True, factor constants from numerator, for example,

        5 * (s**2 + s + 1) / (s**2 + 4)

        See also general, partfrac, standard, timeconst, and ZPK

        """

        var = self.var
        Apoly = self.Apoly
        Bpoly = self.Bpoly
        delay = self.delay

        if factor_const:
            K = sym.cancel(Bpoly.LC() / Apoly.LC())
            if delay != 0:
                K *= sym.exp(self.var * delay)

            # Divide by leading coefficient
            N = Bpoly.monic().as_expr()
            D = Apoly.monic().as_expr()

            if D == 1:
                expr = N
            else:
                expr = sym.Mul(N, 1 / D, evaluate=False)

            if K != 1:
                expr = sym.Mul(K, expr, evaluate=False)
            expr *= self.undef
        else:
            C = Apoly.LC()
            D = Apoly.monic().as_expr()
            N = (Bpoly.as_expr() / C).simplify()
            if D == 1:
                expr = N
            else:
                if N == 1:
                    expr = 1 / D
                else:
                    expr = sym.Mul(N, 1 / D, evaluate=False)
            if delay != 0:
                expr *= sym.exp(self.var * delay)
            expr *= self.undef

        return expr

    def general(self):
        """Convert rational function to general form.

        See also canonical, partfrac, standard, timeconst, and ZPK"""

        B, A, delay, undef = self.as_B_A_delay_undef()

        expr = sym.cancel(B / A, self.var)
        if delay != 0:
            expr *= sym.exp(self.var * delay)

        return expr * undef

    def expandcanonical(self):
        """Expand in terms for different powers with each term
        expressed in canonical form."""

        B, A, delay, undef = self.as_B_A_delay_undef()

        Bpoly = self.Bpoly

        expr = Zero

        for m, c in enumerate(reversed(Bpoly.all_coeffs())):
            term = sym.Mul(c.simplify() * self.var ** m, 1 / A)
            expr += term

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return expr * undef

    def partfrac(self, combine_conjugates=False, damping=None, method=None):
        """Convert rational function into partial fraction form.

        If combine_conjugates is True then the pair of partial
        fractions for complex conjugate poles are combined.

        See also canonical, standard, general, timeconst, and ZPK

        """

        Q, R, F, delay, undef = self.as_QRF(
            combine_conjugates, damping, method)

        result = Q
        for R1, F1 in zip(R, F):
            result += R1 / F1

        if delay != 0:
            result *= sym.exp(-self.var * delay)

        result *= undef
        return result

    def standard(self):
        """Convert rational function into mixed fraction form.

        This is the sum of strictly proper rational function and a
        polynomial.

        See also canonical, general, partfrac, timeconst, and ZPK"""

        Q, M, A, delay, undef = self.as_QMA()

        expr = Q + sym.cancel(M / A, self.var)

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return expr * undef

    def timeconst(self):
        """Convert rational function to time constant form with unity
        lowest power of denominator.

        See also canonical, general, partfrac, standard, and ZPK"""

        B, A, delay, undef = self.as_B_A_delay_undef()

        var = self.var
        Apoly = self.Apoly

        K = Apoly.EC()

        A = A / K
        B = B / K
        return sym.Mul(B, sym.Pow(A, -1), evaluate=False) * sym.exp(self.var * delay) * undef

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
                num *= (var**2 - zeros[0] * var - zeros[1]
                        * var + zeros[0] * zeros[1]).simplify()

        den = 1
        for poles, order in pole_pairs.items():
            for m in range(order):
                den *= (var**2 - poles[0] * var - poles[1]
                        * var + poles[0] * poles[1]).simplify()

        result1 *= (num / den)

        result2 = _zp2tf(zero_singles, pole_singles, 1, var) * undef
        result = K * result1 * result2

        return result

    def residues(self, combine_conjugates=False, damping=None):
        """Return residues of partial fraction expansion.

        This is not much use without the corresponding poles.
        It is better to use as_QRF."""

        Q, R, F, delay, undef = self.as_QRF(combine_conjugates, damping)
        return R

    def coeffs(self):

        Bpoly = self.Bpoly
        Apoly = self.Apoly

        return Bpoly.all_coeffs(), Apoly.all_coeffs()

    @property
    def degree(self):
        """Return the degree (order) of the rational function.

        This the maximum of the numerator and denominator degrees.
        Note zero has a degree of -inf."""

        return max(self.Bpoly.degree(), self.Apoly.degree())

    @property
    def Ndegree(self):
        """Return the degree (order) of the numerator of a rational function.
        Note zero has a degree of -inf."""

        return self.Bpoly.degree()

    @property
    def Ddegree(self):
        """Return the degree (order) of the denominator of a rational function.
        Note zero has a degree of -inf."""

        return self.Apoly.degree()

    @property
    def is_strictly_proper(self):
        """Return True if the degree of the dominator is greater
        than the degree of the numerator."""

        return self.Ddegree > self.Ndegree

    def as_ZPK(self):
        """Decompose expression into zeros, poles, gain, undef where

        expression = K * (prod_n (var - z_n) / (prod_n (var - p_n)) * undef
        """

        Bpoly = self.Bpoly
        Apoly = self.Apoly

        K = sym.cancel(Bpoly.LC() / Apoly.LC())
        if self.delay != 0:
            K *= sym.exp(self.var * self.delay)

        zeros = sym.roots(Bpoly)
        poles = sym.roots(Apoly)

        return zeros, poles, K, self.undef

    def as_QMA(self):
        """Decompose expression into Q, M, A, delay, undef where

        `expression = (Q + M / A) * exp(-delay * var) * undef`"""

        B, A, delay, undef = self.as_B_A_delay_undef()

        # Perform polynomial long division so expr = Q + M / A
        Q, M = sym.div(B, A, self.var)

        return Q, M, A, delay, undef

    def as_QRF_old(self, combine_conjugates=False, damping=None):
        """Decompose expression into Q, R, F, delay, undef where

        expression = (Q + sum_n r_n / f_n) * exp(-delay * var) * undef"""

        Q, M, F, delay, undef = self.as_QMA()

        expr = M / F
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
                        r = sym.simplify(r * (var - pc) **
                                         n + rc * (var - p) ** n)
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

    def _find_residues_ec(self, poles, B):

        from .matrix import matrix_inverse

        # Find residues by solving system of equations (equating coefficients
        # method).  For example,  consider an expression with repeated poles:
        # E(s) = B(S) / (2 * (s - p1)**2 * (s - p2))
        #
        # This can be expressed as partial fractions:
        # E(s) = R_1 / (s - p1)**2 + R_2 / (s - p1) + R_3 / (s - p2)
        # where R_n are the residues we would like to find.
        #
        # The factored denominator expression is
        # Af(s) = ((s - p1)**2 * (s - p2)) = A(s) / C where C is a constant
        # and so the numerator expression is
        # B(s) = E(s) * A(s) = E(s) * Af(s) * C
        #
        # Expanding the partial fractions to have a common denominator gives:
        # B(s) = R_1 * (s - p2) + R_2 * (s - p1) * (s - p2) + R * (s - p1)**2
        #
        # Equating B(s) with E(s) * A(s) and equating powers of s gives
        # system of equation to find the residues.

        var = self.var

        U = []
        F = []
        O = []
        P = []
        i = 1
        A_factored = One
        for pole in poles:
            p = pole.expr
            f = var - p

            for m in range(pole.n):
                u = sym.Symbol('r_%d' % i)
                o = pole.n - m
                U.append(u)
                F.append(f)
                P.append(p)
                O.append(o)
                i += 1
                A_factored *= f

        # This is slow due to use of cancel
        # rhs = Zero
        # for f, o, u in zip(F, O, U):
             # rhs += u * (A_factored / f**o).cancel()

        rhs = Zero
        for i in range(len(U)):
            x = U[i]
            for j in range(len(U)):
                if i == j:
                    continue
                if F[i] is not F[j] or O[j] > O[i]:
                    x *= F[j]
            rhs += x

        lhs = B

        rhspoly = sym.poly(rhs, var)
        lhspoly = sym.poly(lhs, var)

        rc = rhspoly.all_coeffs()
        lc = lhspoly.all_coeffs()
        if len(lc) < len(rc):
            lc = [0] * (len(rc) - len(lc)) + lc

        A, _ = sym.linear_eq_to_matrix(rc, U)

        # Solve system of equations to find residues.
        R = list(matrix_inverse(A) * sym.Matrix(lc))

        return R, P, O

    def _find_residues_sub(self, poles, B):

        var = self.var

        F = []
        O = []
        R = []
        P = []
        M = []

        for pole in poles:
            p = pole.expr
            f = var - p

            for m in range(pole.n):
                o = pole.n - m
                F.append(f)
                P.append(p)
                O.append(o)
                M.append(pole.n)

        for i in range(len(P)):

            if M[i] == O[i]:
                denom = One
                for j in range(len(P)):

                    if i == j:
                        continue
                    if F[i] is not F[j] or O[j] > O[i]:
                        denom *= F[j]
                expr = B / denom
                r = expr.subs(var, P[i])
                R.append(r)
            else:
                expr = expr.diff(var)
                r = expr.subs(var, P[i]) / sym.factorial(M[i] - O[i])
                R.append(r)

        return R, P, O

    def as_QRPO(self, damping=None, method=None):
        """Decompose expression into Q, R, P, O, delay, undef where

        `expression = (Q + sum_n r_n / (var - p_n)**o_n) * exp(-delay * var) * undef`

        `method` can be 'sub' (substitution method, the default) or
        'ec' (equating cofficients method).

        """

        Q, M, A, delay, undef = self.as_QMA()

        expr = M / A
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

        B = sexpr.B
        B /= sexpr.Apoly().LC()

        if method in ('sub', None):
            R, P, O = self._find_residues_sub(poles, B)
        elif method == 'ec':
            R, P, O = self._find_residues_ec(poles, B)
        else:
            raise ValueError('Unknown method ' + method)

        # Remove elements where the residue is zero.
        R, P, O = self._prune_zero_residues(R, P, O)

        return Q, R, P, O, delay, undef

    def _combine_conjugates(self, R, F, P, O):

        Rnew = []
        Fnew = []

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
                Rnew.append(
                    ((self.var - pc) * r + (self.var - p) * rc).expand())
                Fnew.append(((self.var - p) * (self.var - pc)).expand())
            else:
                Rnew.append(r)
                Fnew.append(F[m])
        return Rnew, Fnew

    def as_QRF(self, combine_conjugates=False, damping=None, method=None):
        """Decompose expression into Q, R, F, delay, undef where

        expression = (Q + sum_n r_n / f_n) * exp(-delay * var) * undef"""

        Q, R, P, O, delay, undef = self.as_QRPO(damping, method)

        F = [(self.var - p)**o for p, o in zip(P, O)]

        if combine_conjugates:
            R, F = self._combine_conjugates(R, F, P, O)

        return Q, R, F, delay, undef

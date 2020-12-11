"""
This module provides support for rational functions.

Copyright 2016--2020 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from sympy.core.mul import _unevaluated_Mul as uMul
from sympy.core.add import _unevaluated_Add as uAdd
from .sym import sympify, AppliedUndef

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
        scale = sym.S.One
    
        for factor in expr.as_ordered_factors():
            if (factor.is_Pow and factor.args[1] == sym.S.One / 2 and
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
    
        offset = sym.S.Zero
        scale = sym.S.One
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

        scale = sym.S.One
        scale2 = sym.S.One    
        offset = sym.S.One
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
            dexpr = sym.S.Zero
        return scale, offset, scale2, dexpr

    @property
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

    N = sym.S.One
    D = sym.S.One    
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
    
    N = sym.S.One
    D = sym.S.One    
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
        
    return uMul(K, *(zz + pp))


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
    return uMul(K, *(zz + pp))


def _pr2tf(poles, residues, var=None):
    """Create a transfer function from lists of poles and residues.

    """
    poles = sympify(poles)
    residues = sympify(residues)

    return uAdd(*[r / (var - p) for r, p in zip(residues, poles)])


def as_ratfun_delay(expr, var):
    delay = sym.S.Zero

    if expr.is_rational_function(var):
        N, D = expr.as_numer_denom()
        return N, D, delay

    F = sym.factor(expr).as_ordered_factors()

    rf = sym.S.One
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
    delay = sym.S.Zero
    undef = sym.S.One
    
    if expr.is_rational_function(var):
        N, D = as_numer_denom(expr, var)
        return N, D, delay, undef

    F = sym.factor(expr).as_ordered_factors()

    rf = sym.S.One
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
        
        const = sym.S.One
        undef = sym.S.One
        rest = sym.S.One        
    
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
    
    def roots(self):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        return sym.roots(sym.Poly(self.expr, self.var))

    def zeros(self):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        return Ratfun(self.numerator, self.var).roots()    

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
        return m2


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

            expr = sym.Mul(K, expr, undef, evaluate=False)
        else:
            C = Dpoly.LC()
            D = Dpoly.monic().as_expr()
            N = (Npoly.as_expr() / C).simplify()
            if D == 1:
                expr = N
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
        
        expr = sym.S.Zero

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

    def ZPK(self):
        """Convert to zero-pole-gain (ZPK) form.

        See also canonical, general, standard, timeconst, and partfrac"""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        var = self.var        
        Npoly = sym.Poly(N, var)
        Dpoly = sym.Poly(D, var)
        
        K = sym.cancel(Npoly.LC() / Dpoly.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        zeros = sym.roots(Npoly)
        poles = sym.roots(Dpoly)

        return _zp2tf(zeros, poles, K, self.var) * undef

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

    def as_QMD(self):
        """Decompose expression into Q, M, D, delay, undef where

        expression = (Q + M / D) * exp(-delay * var) * undef"""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        # Perform polynomial long division so expr = Q + M / D
        Q, M = sym.div(N, D, self.var)

        return Q, M, D, delay, undef
        
    def as_QRD(self, combine_conjugates=False, damping=None):
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
                        r = sym.limit(
                            sym.diff(expr2, var, m), var, p) / sym.factorial(m)
                        rc = r.conjugate()
                        r = sym.simplify(r * (var - pc) ** n + rc * (var - p) ** n)
                        R.append(r)
                        D.append(D2 ** n)
            else:
                D2 = var - p

                if o == 1:
                    r = sexpr.residue(p, poles)
                    R.append(r)
                    D.append(D2)
                else:
                    # Handle repeated real poles.
                    expr2 = expr * (var - p) ** o
                    for n in range(1, o + 1):
                        m = o - n
                        r = sym.limit(
                            sym.diff(expr2, var, m), var, p) / sym.factorial(m)

                        R.append(r)
                        D.append(D2 ** n)                        
                                   
        return Q, R, D, delay, undef
            

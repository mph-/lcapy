"""
This module provides support for rational functions.

Copyright 2016--2019 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from sympy.core.mul import _unevaluated_Mul as uMul
from .sym import sympify


def _zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists of zeros and poles,
    and from a constant gain"""

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
        N, D = expr.as_numer_denom()
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
        if isinstance(f, sym.function.AppliedUndef):
            undef *= f
            continue
                
        rf *= f

    if not rf.is_rational_function(var):
        raise ValueError('Expression not a product of rational function'
                         ' exponential, and undefined functions')
    
    N, D = rf.as_numer_denom()
    return N, D, delay, undef


def as_residue_parts(expr, var):
        
    N, D, delay = as_ratfun_delay(expr, var)

    # Perform polynomial long division so expr = Q + M / D
    Q, M = sym.div(N, D, var)
    expr = M / D
        
    sexpr = Ratfun(expr, var)

    P = sexpr.poles()
    F = []
    R = []
    for p in P:

        # Number of occurrences of the pole.
        N = P[p]

        f = var - p

        if N == 1:
            F.append(f)
            R.append(sexpr.residue(p, P))
            continue

        # Handle repeated poles.
        expr2 = expr * f ** N
        for n in range(1, N + 1):
            m = N - n
            F.append(f ** n)
            dexpr = sym.diff(expr2, var, m)
            R.append(sym.limit(dexpr, var, p) / sym.factorial(m))

    return F, R, Q, delay


class Ratfun(object):

    def __init__(self, expr, var):
        self.expr = expr
        self.var = var

    def as_residue_parts(self):
        """Return residues of expression"""

        return as_residue_parts(self.expr, self.var)

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
            if isinstance(f, sym.function.AppliedUndef):
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

    def poles(self):
        """Return poles of expression as a dictionary
        Note this may not find them all."""

        return Ratfun(self.denominator, self.var).roots()    

    def residue(self, pole, poles):

        expr = self.expr
        var = self.var
        
        # Remove pole from list of poles; sym.cancel
        # doesn't always work, for example, for complex poles.
        poles2 = poles.copy()
        poles2[pole] -= 1
        
        numer, denom = expr.as_numer_denom()
        Dpoly = sym.Poly(denom, var)
        K = Dpoly.LC()
        
        D = [(var - p) ** poles2[p] for p in poles2]
        denom = sym.Mul(K, *D)
        
        d = sym.limit(denom, var, pole)
        
        if d != 0:
            tmp = numer / denom
            return sym.limit(tmp, var, pole)

        # Use l'Hopital's rule
        tmp = numer / denom
        tmp = sym.diff(tmp, var)
        
        return sym.limit(tmp, var, pole)

    @property
    def numerator(self):
        """Return numerator of rational function"""

        numer, denom = self.expr.as_numer_denom()
        return numer

    @property
    def denominator(self):
        """Return denominator of rational function"""

        numer, denom = self.expr.as_numer_denom()
        return denom

    def canonical(self):
        """Convert rational function to canonical form with unity
        highest power of denominator.

        See also general, partfrac, mixedfrac, and ZPK"""

        try:
            N, D, delay, undef = self.as_ratfun_delay_undef()
        except ValueError:
            # TODO: copy?
            return self.expr

        var = self.var        
        Dpoly = sym.Poly(D, var)
        Npoly = sym.Poly(N, var)

        K = sym.cancel(Npoly.LC() / Dpoly.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        # Divide by leading coefficient
        Nm = Npoly.monic()
        Dm = Dpoly.monic()

        expr = K * (Nm / Dm)

        return expr * undef

    def general(self):
        """Convert rational function to general form.

        See also canonical, partfrac, mixedfrac, and ZPK"""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        expr = sym.cancel(N / D, self.var)
        if delay != 0:
            expr *= sym.exp(self.var * delay)

        return expr * undef

    def partfrac(self):
        """Convert rational function into partial fraction form.

        See also canonical, mixedfrac, general, and ZPK"""

        N, D, delay, undef = self.as_ratfun_delay_undef()
        expr = N / D
        
        try:
            F, R, Q, delay2 = as_residue_parts(expr, self.var)
            
        except ValueError:
            # Try splitting into terms
            expr = expr.expand()
            if not expr.is_Add:
                raise ValueError('Cannot convert to partial fraction')
            result = 0
            for arg in expr.args:
                result += Ratfun(arg, self.var).partfrac()
            return result * undef

        expr = Q.as_expr()
        for f, r in zip(F, R):
            expr += r / f

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return expr * undef

    def mixedfrac(self):
        """Convert rational function into mixed fraction form.

        See also canonical, general, partfrac and ZPK"""

        N, D, delay, undef = self.as_ratfun_delay_undef()
        var = self.var        

        # Perform polynomial long division so expr = Q + M / D        
        Q, M = sym.div(N, D, var)
        expr = Q + sym.cancel(M / D, var)

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return expr * undef

    def timeconst(self):
        """Convert rational function to time constant form with unity
        lowest power of denominator.

        See also canonical, general, partfrac, mixedfrac, and ZPK"""

        try:
            N, D, delay, undef = self.as_ratfun_delay_undef()
        except ValueError:
            # TODO: copy?
            return self.expr

        K = sym.cancel(N.EC() / D.EC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        # Divide by leading coefficient
        Nm = (N / N.EC()).simplify()
        Dm = (D / D.EC()).simplify()
        
        expr = K * (Nm / Dm)

        return expr * undef

    def ZPK(self):
        """Convert to pole-zero-gain (PZK) form.

        See also canonical, general, mixedfrac, and partfrac"""

        N, D, delay, undef = self.as_ratfun_delay_undef()

        var = self.var        
        Dpoly = sym.Poly(D, var)
        Npoly = sym.Poly(N, var)
        
        K = sym.cancel(Npoly.LC() / Dpoly.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        zeros = sym.roots(Npoly)
        poles = sym.roots(Dpoly)

        return _zp2tf(zeros, poles, K, self.var) * undef

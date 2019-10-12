from .expr import Expr, expr
import sympy as sym

class sfwExpr(Expr):

    def __init__(self, val, **assumptions):

        super(sfwExpr, self).__init__(val, **assumptions)

    def roots(self, aslist=False):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        roots = self._ratfun.roots()
        if not aslist:
            return expr(roots)
        rootslist = []
        for root, count in roots.items():
            rootslist += [root] * count
        return expr(rootslist)
            
    def zeros(self, aslist=False):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        return self.N.roots(aslist)

    def poles(self, aslist=False):
        """Return poles of expression as a dictionary
        Note this may not find them all."""

        return self.D.roots(aslist)

    def residues(self):
        """Return list of residues of partial fraction expansion."""

        return expr(self._ratfun.residues())

    def canonical(self, factor_const=True):
        """Convert rational function to canonical form (standard form) with
        unity highest power of denominator.  For example,

        5 * (s**2 + s + 1) / (s**2 + 4)

        See also general, partfrac, mixedfrac, timeconst, and ZPK

        """

        return self.__class__(self._ratfun.canonical(factor_const), **self.assumptions)

    def general(self):
        """Convert rational function to general form.  For example,

        (5 * s**2 + 10 * s + 5) / (s**2 + 4)

        See also canonical, partfrac, mixedfrac, timeconst, and ZPK."""

        return self.__class__(self._ratfun.general(), **self.assumptions)

    def partfrac(self):
        """Convert rational function into partial fraction form.   For example,

        5 + (5 - 15 * j / 4) / (s + 2 * j) + (5 + 15 * j / 4) / (s - 2 * j)

        See also canonical, mixedfrac, general, timeconst, and ZPK."""

        return self.__class__(self._ratfun.partfrac(), **self.assumptions)

    def mixedfrac(self):
        """Convert rational function into mixed fraction form.  For example,

        (5 * s - 5) / (s**2 + 4) + 5

        See also canonical, general, partfrac, timeconst, and ZPK."""

        return self.__class__(self._ratfun.mixedfrac(), **self.assumptions)

    def timeconst(self):
        """Convert rational function into time constant form.  For example,

        5 * (s**2 + 2 * s + 1) / (4 * (s**2 / 4 + 1))

        See also canonical, general, mixedfrac, partfrac and ZPK."""

        return self.__class__(self._ratfun.timeconst(), **self.assumptions)

    def ZPK(self):
        """Convert to zero-pole-gain (ZPK) form (factored form).  For example,

        5 * (s + 1)**2 / ((s - 2 * j) * (s + 2 * j))

        See also canonical, general, mixedfrac, partfrac, and timeconst.

        """

        return self.__class__(self._ratfun.ZPK(), **self.assumptions)

    def expandcanonical(self):
        """Expand in terms for different powers with each term
        expressed in canonical form.  For example,

        s / (s**2 + 4) + 5 / (s**2 + 4)

        See also canonical, general, partfrac, timeconst, and ZPK."""

        return self.__class__(self._ratfun.expandcanonical(), **self.assumptions)

    def coeffs(self):
        """Return list of coeffs assuming the expr is a polynomial in s.  The
        highest powers come first.  This will fail for a rational function.
        Instead use expr.N.coeffs or expr.D.coeffs for numerator
        or denominator respectively."""

        try:
            z = sym.Poly(self.expr, self.var)
        except:
            raise ValueError('Use .N or .D attribute to specify numerator or denominator of rational function')
        return expr(z.all_coeffs())

    def normcoeffs(self):
        """Return list of coeffs (normalised so the highest power is 1)
        assuming the expr is a polynomial in s.  The highest powers
        come first.  This will fail for a rational function.  Instead
        use expr.N.normcoeffs or expr.D.normcoeffs for numerator or
        denominator respectively."""

        try:
            z = sym.Poly(self.expr, self.var)
        except:
            raise ValueError('Use .N or .D attribute to specify numerator or denominator of rational function')            
        c = z.all_coeffs()
        return expr([sym.simplify(c1 / c[0]) for c1 in c])

from .expr import Expr, ExprDict, ExprList
import sympy as sym

class sfwExpr(Expr):

    def __init__(self, val, **assumptions):

        super(sfwExpr, self).__init__(val, **assumptions)

    def roots(self):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        return ExprDict(self._ratfun.roots())

    def zeros(self):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        return self.N.roots()

    def poles(self):
        """Return poles of expression as a dictionary
        Note this may not find them all."""

        return self.D.roots()

    def residues(self):
        """Return list of residues of partial fraction expansion."""

        return ExprList(self._ratfun.residues())

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

        z = sym.Poly(self.expr, self.var)
        return z.all_coeffs()

    def normcoeffs(self):
        """Return list of coeffs (normalised so the highest power is 1)
        assuming the expr is a polynomial in s.  The highest powers
        come first.  This will fail for a rational function.  Instead
        use expr.N.normcoeffs or expr.D.normcoeffs for numerator or
        denominator respectively."""
        
        z = sym.Poly(self.expr, self.var)
        c = z.all_coeffs()
        return [sym.simplify(c1 / c[0]) for c1 in c]

from .expr import Expr, ExprDict

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

    def canonical(self, factor_const=True):
        """Convert rational function to canonical form with unity
        highest power of denominator.   For example,

        5 * (s**2 + s + 1) / (s**2 + 4)

        See also general, partfrac, mixedfrac, timeconst, and ZPK"""

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
        """Convert to zero-pole-gain (ZPK) form.  For example,

        5 * (s + 1)**2 / ((s - 2 * j) * (s + 2 * j))

        See also canonical, general, mixedfrac, partfrac, and timeconst."""

        return self.__class__(self._ratfun.ZPK(), **self.assumptions)

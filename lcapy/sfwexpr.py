from .expr import Expr, Exprdict

class sfwExpr(Expr):

    def __init__(self, val, **assumptions):

        super(sfwExpr, self).__init__(val, **assumptions)

    def roots(self):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        return Exprdict(self._ratfun.roots())

    def zeros(self):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        return self.N.roots()

    def poles(self):
        """Return poles of expression as a dictionary
        Note this may not find them all."""

        return self.D.roots()

    def canonical(self):
        """Convert rational function to canonical form with unity
        highest power of denominator.

        See also general, partfrac, mixedfrac, timeconst, and ZPK"""

        return self.__class__(self._ratfun.canonical(), **self.assumptions)

    def general(self):
        """Convert rational function to general form.

        See also canonical, partfrac, mixedfrac, timeconst, and ZPK."""

        return self.__class__(self._ratfun.general(), **self.assumptions)

    def partfrac(self):
        """Convert rational function into partial fraction form.

        See also canonical, mixedfrac, general, timeconst, and ZPK."""

        return self.__class__(self._ratfun.partfrac(), **self.assumptions)

    def mixedfrac(self):
        """Convert rational function into mixed fraction form.

        See also canonical, general, partfrac, timeconst, and ZPK."""

        return self.__class__(self._ratfun.mixedfrac(), **self.assumptions)

    def timeconst(self):
        """Convert rational function into time constant form.

        See also canonical, general, mixedfrac, partfrac and ZPK."""

        return self.__class__(self._ratfun.timeconst(), **self.assumptions)

    
    def ZPK(self):
        """Convert to pole-zero-gain (PZK) form.

        See also canonical, general, mixedfrac, partfrac, and timeconst."""

        return self.__class__(self._ratfun.ZPK(), **self.assumptions)

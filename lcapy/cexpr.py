from .expr import Expr
from .sym import symbols_find

__all__ = ('Vconst', 'Iconst')


class cExpr(Expr):

    """Constant real expression or symbol.

    If symbols in the expression are known to be negative, use
    cExpr(expr, positive=False)

    """

    def __init__(self, val, **assumptions):

        symbols = symbols_find(val)
        for symbol in ('s', 'omega', 't', 'f'):
            if symbol in symbols:
                raise ValueError(
                    'constant expression %s cannot depend on %s' % (val, symbol))

        assumptions['dc'] = True        
        super(cExpr, self).__init__(val, **assumptions)

    def rms(self):
        return {Vconst: Vt, Iconst : It}[self.__class__](self)

    def laplace(self):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    def canonical(self, factor_const=True):
        # Minor optimisation
        return self

    def time(self, **assumptions):
        """Convert to time domain."""
        
        return tExpr(self)

class Vconst(cExpr):

    superkind = 'Voltage'
    
    def __init__(self, val, **assumptions):

        super(Vconst, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Vt

    def cpt(self):
        from .oneport import Vdc
        return Vdc(self)

    def time(self, **assumptions):
        return Vt(self)

    
class Iconst(cExpr):

    superkind = 'Current'
    
    def __init__(self, val, **assumptions):

        super(Iconst, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = It

    def cpt(self):
        from .oneport import Idc
        return Idc(self)

    def time(self, **assumptions):
        return It(self)


from .texpr import t, It, Vt, tExpr
from .sexpr import s

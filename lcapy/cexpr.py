from .expr import Expr
from .sympify import symbols_find

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

        super(cExpr, self).__init__(val, **assumptions)

    def rms(self):
        return {Vconst: Vt, Iconst : It}[self.__class__](self)

    def laplace(self):
        """Convert to Laplace domain representation."""

        return self.time().laplace()


class Vconst(cExpr):

    def __init__(self, val, **assumptions):

        assumptions['dc'] = True        
        super(Vconst, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Vt

    def cpt(self):
        from .oneport import Vdc
        return Vdc(self)

    def time(self, **assumptions):
        return Vt(self)

    
class Iconst(cExpr):

    def __init__(self, val, **assumptions):

        super(Iconst, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = It

    def cpt(self):
        from .oneport import Idc
        return Idc(self)

    def time(self, **assumptions):
        return It(self)


from .texpr import Ht, It, Vt, Yt, Zt, tExpr    


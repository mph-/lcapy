"""This module implements the sMatrix class for a matrix of
Laplace-domain expressions.

Copyright 2019 Michael Hayes, UCECE

"""

from .matrix import Matrix
from .laplace import inverse_laplace_transform


class sMatrix(Matrix):
    from .sexpr import sExpr    
    _typewrap = sExpr

    def inverse_laplace(self, **assumptions):

        def ilt(expr):
            from .sym import ssym, tsym
            return inverse_laplace_transform(expr, ssym, tsym, **assumptions)
        
        return tMatrix(self.applyfunc(ilt))

    def canonical(self):

        return self.applyfunc(self._typewrap.canonical)

    def general(self):

        return self.applyfunc(self._typewrap.general)

    def mixedfrac(self):

        return self.applyfunc(self._typewrap.mixedfrac)

    def partfrac(self):

        return self.applyfunc(self._typewrap.partfrac)

    def timeconst(self):

        return self.applyfunc(self._typewrap.timeconst)   

    def ZPK(self):

        return self.applyfunc(self._typewrap.ZPK)    
    

from .tmatrix import tMatrix
    

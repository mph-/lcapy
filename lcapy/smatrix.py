from .matrix import Matrix
from .laplace import inverse_laplace_transform


class sMatrix(Matrix):
    from .sexpr import sExpr    
    _typewrap = sExpr

    def inverse_laplace(self, **assumptions):

        def ilt(expr):
            from .sym import ssym, tsym
            return inverse_laplace_transform(expr, ssym, tsym, **assumptions)
        
        return self.applyfunc(ilt)    

    

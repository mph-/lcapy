"""This module implements the TimeDomainMatrix class for a matrix of
time-domain expressions.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from .matrix import Matrix
from .laplace import laplace_transform

class TimeDomainMatrix(Matrix):
    from .texpr import TimeDomainExpression    
    _typewrap = TimeDomainExpression    

    def laplace(self):

        def lt(expr):
            from .sym import ssym, tsym
            return laplace_transform(expr, tsym, ssym)
        
        return LaplaceDomainMatrix(self.applyfunc(lt))
    
from .smatrix import LaplaceDomainMatrix

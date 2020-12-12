"""This module implements the LaplaceDomainMatrix class for a matrix of
Laplace-domain expressions.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from .matrix import Matrix
from .laplace import inverse_laplace_transform


class LaplaceDomainMatrix(Matrix):
    from .sexpr import LaplaceDomainExpression    
    _typewrap = LaplaceDomainExpression

    def inverse_laplace(self, **assumptions):

        def ilt(expr):
            from .sym import ssym, tsym
            return inverse_laplace_transform(expr, ssym, tsym, **assumptions)
        
        return TimeDomainMatrix(self.applyfunc(ilt))

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
    

class LaplaceDomainVoltageMatrix(LaplaceDomainMatrix):
    from .sexpr import LaplaceDomainVoltage    
    _typewrap = LaplaceDomainVoltage

    
class LaplaceDomainCurrentMatrix(LaplaceDomainMatrix):
    from .sexpr import LaplaceDomainCurrent    
    _typewrap = LaplaceDomainCurrent

    
class LaplaceDomainAdmittanceMatrix(LaplaceDomainMatrix):
    from .sexpr import LaplaceDomainAdmittance    
    _typewrap = LaplaceDomainAdmittance

    
class LaplaceDomainImpedanceMatrix(LaplaceDomainMatrix):
    from .sexpr import LaplaceDomainImpedance    
    _typewrap = LaplaceDomainImpedance    
    
    
from .tmatrix import TimeDomainMatrix
    


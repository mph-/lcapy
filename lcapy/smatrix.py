"""This module implements the LaplaceDomainMatrix class for a matrix of
Laplace-domain expressions.

Copyright 2019--2023 Michael Hayes, UCECE

"""

from .matrix import Matrix


class LaplaceDomainMatrix(Matrix):
    from .exprclasses import LaplaceDomainExpression
    _typewrap = LaplaceDomainExpression

    def ILT(self, **assumptions):

        def func(expr):
            return expr.ILT(**assumptions)

        return TimeDomainMatrix(self.applyfunc(func))


class LaplaceDomainVoltageMatrix(LaplaceDomainMatrix):
    from .exprclasses import LaplaceDomainVoltage
    _typewrap = LaplaceDomainVoltage


class LaplaceDomainCurrentMatrix(LaplaceDomainMatrix):
    from .exprclasses import LaplaceDomainCurrent
    _typewrap = LaplaceDomainCurrent


class LaplaceDomainAdmittanceMatrix(LaplaceDomainMatrix):
    from .exprclasses import LaplaceDomainAdmittance
    _typewrap = LaplaceDomainAdmittance


class LaplaceDomainImpedanceMatrix(LaplaceDomainMatrix):
    from .exprclasses import LaplaceDomainImpedance
    _typewrap = LaplaceDomainImpedance


from .tmatrix import TimeDomainMatrix  # nopep8

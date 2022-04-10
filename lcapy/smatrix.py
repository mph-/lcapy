"""This module implements the LaplaceDomainMatrix class for a matrix of
Laplace-domain expressions.

Copyright 2019--2021 Michael Hayes, UCECE

"""

from .matrix import Matrix


class LaplaceDomainMatrix(Matrix):
    from .sexpr import LaplaceDomainExpression
    _typewrap = LaplaceDomainExpression

    def ILT(self, **assumptions):

        def func(expr):
            return expr.ILT(**assumptions)

        return TimeDomainMatrix(self.applyfunc(func))


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


from .tmatrix import TimeDomainMatrix  # nopep8

"""This module implements the DiscreteTimeDomainMatrix class for a matrix of
time-domain expressions.

Copyright 2021 Michael Hayes, UCECE

"""

from .matrix import Matrix


class DiscreteTimeDomainMatrix(Matrix):
    from .nexpr import DiscreteTimeDomainExpression
    _typewrap = DiscreteTimeDomainExpression

    def ZT(self):

        def func(expr):
            return expr.ZT()

        return ZDomainMatrix(self.applyfunc(func))


from .zmatrix import ZDomainMatrix  # nopep8

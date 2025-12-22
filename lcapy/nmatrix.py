"""This module implements the DiscreteTimeDomainMatrix class for a matrix of
time-domain expressions.

Copyright 2021--2025 Michael Hayes, UCECE

"""

from .matrix import Matrix


class DiscreteTimeDomainMatrix(Matrix):
    from .nexpr import DiscreteTimeDomainExpression
    _typewrap = DiscreteTimeDomainExpression

    def ZT(self):

        def func(expr):
            return DiscreteTimeDomainExpression(expr).ZT()

        return ZDomainMatrix(Matrix(self).applyfunc(func))


from .zmatrix import ZDomainMatrix  # nopep8
from .nexpr import DiscreteTimeDomainExpression # nopep

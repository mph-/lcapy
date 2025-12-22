"""This module implements the TimeDomainMatrix class for a matrix of
time-domain expressions.

Copyright 2019--2025 Michael Hayes, UCECE

"""

from .matrix import Matrix


class TimeDomainMatrix(Matrix):
    from .texpr import TimeDomainExpression
    _typewrap = TimeDomainExpression

    def LT(self):

        def func(expr):
            return TimeDomainExpression(expr).LT()

        return LaplaceDomainMatrix(Matrix(self).applyfunc(func))


from .smatrix import LaplaceDomainMatrix  # nopep8
from .texpr import TimeDomainExpression # nopep

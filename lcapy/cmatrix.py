"""This module implements the ConstantDomainMatrix class for a matrix of
time-domain expressions.

Copyright 2023--2025 Michael Hayes, UCECE

"""

from .matrix import Matrix


class ConstantDomainMatrix(Matrix):
    from .cexpr import ConstantDomainExpression
    _typewrap = ConstantDomainExpression

"""This module provides inductance support.

Copyright 2021--2022 Michael Hayes, UCECE

"""
from .cexpr import cexpr
from .units import u as uu
from warnings import warn


def inductance(arg, **assumptions):

    expr1 = cexpr(arg, frequency=True, **assumptions)
    if expr1.is_imaginary:
        warn('Inductance %s should be real' % expr1)

    expr1.units = uu.henry
    return expr1

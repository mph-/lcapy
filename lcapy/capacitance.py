"""This module provides capacitance support.

Copyright 2021 Michael Hayes, UCECE

"""
from .cexpr import cexpr
from .units import u as uu
from warnings import warn


def capacitance(arg, **assumptions):

    expr1 = cexpr(arg, frequency=True, **assumptions)
    if expr1.is_imaginary:
        warn('Capacitance %s should be real' % expr1)

    expr1.units = uu.farad
    return expr1

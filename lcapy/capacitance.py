"""This module provides capacitance support.

Copyright 2021 Michael Hayes, UCECE

"""
from .cexpr import cexpr
from .units import u as uu

def capacitance(arg, **assumptions):

    expr1 = cexpr(arg, **assumptions)
    expr1.units = uu.farad

    return expr1

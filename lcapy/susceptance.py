"""This module provides susceptance support.

Copyright 2021 Michael Hayes, UCECE

"""
from .cexpr import cexpr
from .units import u as uu

def susceptance(arg, **assumptions):

    expr1 = cexpr(arg, **assumptions)
    expr1.units = uu.siemens

    return expr1

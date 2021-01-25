"""This module provides resistance support.

Copyright 2021 Michael Hayes, UCECE

"""
from .cexpr import cexpr
from .units import u as uu

def resistance(arg, **assumptions):

    # Perhaps, relax to expr?
    expr1 = cexpr(arg, **assumptions)
    expr1.units = uu.ohm

    return expr1

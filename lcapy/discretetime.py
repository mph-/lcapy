"""This module provides discrete-time support.

It introduces three special variables:
   n for discrete-time sequences
   k for discrete-frequency sequences
   z for z-transforms.

Copyright 2020 Michael Hayes, UCECE

"""

import sympy as sym
from .sym import sympify
from .nexpr import nexpr, nExpr, n
from .kexpr import kexpr, kExpr, k
from .zexpr import zexpr, zExpr, z
from .dsym import nsym, ksym, zsym, dt, df

from .expr import expr as expr1
from .transform import transform as transform1
from .transform import call as call1
from .functions import Function
from .ztransform import *
from .seq import seq


def expr(arg, **assumptions):

    expr = sympify(arg, **assumptions)

    symbols = expr.free_symbols    

    if nsym in symbols:
        return nexpr(expr, **assumptions)
    elif ksym in symbols:
        return kexpr(expr, **assumptions)    
    elif zsym in symbols:
        return zexpr(expr, **assumptions)    

    return expr1(arg, **assumptions)


def transform(expr, arg, **assumptions):

    # Is this wise?   It makes sense for Voltage and Impedance objects
    # but may cause too much confusion for other expressions
    if arg is n and isinstance(expr, zExpr):
        return expr.IZT(**assumptions)
    elif arg is n and isinstance(expr, kExpr):
        return expr.IDFT(**assumptions)        
    elif arg is z and isinstance(expr, nExpr):
        return expr.ZT(**assumptions)
    elif arg is z and isinstance(expr, kExpr):
        return expr.IDFT(**assumptions).ZT(**assumptions)
    elif arg is k and isinstance(expr, nExpr):
        return expr.DFT(**assumptions)
    elif arg is k and isinstance(expr, zExpr):
        return expr.IZT(**assumptions).DFT(**assumptions)

    return transform1(expr, arg, **assumptions)    


def call(expr, arg, **assumptions):

    if id(arg) in (id(n), id(z), id(k)):
        return transform(expr, arg, **assumptions)

    if arg in (n, k, z):
        return transform(expr, arg, **assumptions)    
    
    return call1(expr, arg, **assumptions)


print('Warning, this is experimental and probably riddled with bugs!')

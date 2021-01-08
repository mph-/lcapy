"""This module provides discrete-time support.

It introduces three special variables:
   n for discrete-time sequences
   k for discrete-frequency sequences
   z for z-transforms.

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy as sym
from .sym import sympify
from .nexpr import nexpr, n
from .kexpr import kexpr, k
from .zexpr import zexpr, z
from .dsym import nsym, ksym, zsym, dt, df

from .expr import expr as expr1
from .transform import transform as transform1
from .transform import call as call1
from .functions import Function
from .ztransform import *
from .seq import seq


def expr(arg, **assumptions):

    # Handle container args.
    if not isinstance(arg, str) and hasattr(arg, '__iter__'):
        return expr1(arg, **assumptions)
    
    expr = sympify(arg, **assumptions)

    symbols = expr.free_symbols    

    if nsym in symbols:
        return nexpr(expr, **assumptions)
    elif ksym in symbols:
        return kexpr(expr, **assumptions)    
    elif zsym in symbols:
        return zexpr(expr, **assumptions)    

    return expr1(arg, **assumptions)


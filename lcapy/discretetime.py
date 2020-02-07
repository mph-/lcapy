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

from lcapy import expr as expr1
from .functions import Function


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

# TODO:
# Add method to create difference equations, say difference_equations(output='y', input='x')
# Approximate tExpr with nExpr
# Rewrite rational functions in terms of z^-1
# Symbol for z^-1, say invz or iz?  Would need special casing to handle invz * z etc.
# Transforms, perhaps have a Transformer class with registered transforms?
# Unit lag operator L() ? 
# Difference operator D() for first-order differences
# Perhaps transforms ZT(), IZT(), FT(), IFT(), LT(), ILT() ?
# Fix summations over unit impulses


print('Warning this is experimental and probably riddled with bugs')

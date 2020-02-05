from .sym import sympify
from .nexpr import nexpr, nExpr, n
from .kexpr import kexpr, kExpr, k
from .zexpr import zexpr, zExpr, z
from .dsym import nsym, ksym, zsym

from lcapy import expr as expr1

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
# add method to create difference equations, say difference_equations(output='y', input='x')
# Approximate tExpr with nExpr
# Add kExpr for discrete frequency domain as conjugate to nExpr
# Impulse function
# Z transforms
# Sequences -> sum of delayed impulses
# Rewrite rational functions in terms of z^-1
# Symbol for z^-1, say invz?  Would need special casing to handle invz * z etc.
# Transforms, perhaps have a Transformer class with registered transforms?



"""This module performs transformations between domains.

Copyright 2018--2020 Michael Hayes, UCECE

"""

from .sym import sympify, pi
from .symbols import f, s, t, omega, jomega

def transform(expr, arg, **assumptions):
    """If arg is f, s, t, omega, jomega perform domain transformation,
    otherwise perform substitution.

    Note (5 * s)(omega) will fail since 5 * s is assumed not to be
    causal and so Fourier transform is unknown.  However, Zs(5 *
    s)(omega) will work since Zs is assumed to be causal.

    """

    # handle things like (3 * s)(5 * s)
    if isinstance(expr, arg.__class__) and not isinstance(expr, Super):
        return expr.subs(arg)

    # Handle expr(t), expr(s), expr(f)
    if arg is t:
        return expr.time(**assumptions)
    elif arg is s:
        return expr.laplace(**assumptions)
    elif arg is f:
        if isinstance(expr, omegaExpr):
            return expr.subs(2 * pi * f)
        else:
            return expr.fourier(**assumptions)    

    # Handle expr(texpr), expr(sexpr), expr(fexpr).  For example,
    # expr(2 * f).
    result = None 
    if isinstance(arg, tExpr):
        result = expr.time(**assumptions)
    elif isinstance(arg, sExpr):
        result = expr.laplace(**assumptions)
    elif isinstance(arg, fExpr):
        if isinstance(expr, omegaExpr):
            result = expr.subs(2 * pi * f)
        else:
            result = expr.fourier(**assumptions)

    if result is not None:
        return result.__class__(result.subs(arg))

    sarg = arg
    try:
        sarg.has(t)
    except:
        sarg = sympify(arg)            

    if sarg.has(jomega):
        if isinstance(expr, sExpr):
            result = expr.subs(arg)
            return result
        else:
            result = expr.laplace(**assumptions)
    elif isinstance(arg, omegaExpr):
        if isinstance(expr, fExpr):
            result = expr.subs(omega / (2 * pi))
        elif isinstance(expr, cExpr):
            result = expr
        else:
            result = expr.fourier(**assumptions).subs(omega / (2 * pi))
    elif sarg.is_constant():
        if not isinstance(expr, Super):
            result = expr.time(**assumptions)
        else:
            result = expr
    else:
        raise ValueError('Can only return t, f, s, or omega domains')
    
    return result.subs(arg)


def call(expr, arg, **assumptions):

    if id(arg) in (id(f), id(s), id(t), id(omega), id(jomega)):
        return transform(expr, arg, **assumptions)

    if arg in (f, s, t, omega, jomega):
        return transform(expr, arg, **assumptions)    
    
    return expr.subs(arg)


from .cexpr import cExpr    
from .fexpr import fExpr    
from .sexpr import sExpr
from .texpr import tExpr
from .omegaexpr import omegaExpr
from .super import Super


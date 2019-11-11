from .sym import sympify, pi
from .symbols import f, s, t, omega, jomega

def transform(expr, arg, **assumptions):

    if isinstance(expr, arg.__class__) and not isinstance(expr, Super):
        return expr.subs(arg)

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
        return result.subs(arg)

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


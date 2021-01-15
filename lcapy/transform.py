"""This module performs transformations between domains.

Copyright 2018--2020 Michael Hayes, UCECE

"""

from .sym import sympify, pi
from .symbols import f, s, t, omega, j, jw, jw0
from .expr import expr as expr1
from .expr import Expr


def transform1(expr, arg, **assumptions):

    # Handle things like (3 * s)(5 * s)
    if isinstance(expr, arg.__class__) and not isinstance(expr, Superposition):
        return expr.subs(arg)

    # Handle expr(t), expr(s), expr(f), expr(omega)
    if arg is t:
        return expr.time(**assumptions)
    elif arg is s:
        return expr.laplace(**assumptions)
    elif arg is f:
        return expr.fourier(**assumptions)
    elif arg is omega:
        return expr.angular_fourier(**assumptions)
    elif arg.has(j):
        return expr.phasor(omega=arg / j, **assumptions)    

    # Handle expr(texpr), expr(sexpr), expr(fexpr), expr(omegaexpr).  For example,
    # expr(2 * f).
    result = None 
    if isinstance(arg, TimeDomainExpression):
        result = expr.time(**assumptions)
    elif isinstance(arg, LaplaceDomainExpression):
        result = expr.laplace(**assumptions)
    elif isinstance(arg, FourierDomainExpression):
        result = expr.fourier(**assumptions)
    elif isinstance(arg, AngularFourierDomainExpression):
        result = expr.angular_fourier(**assumptions)        
    elif arg.has(j):
        result = expr.phasor(omega=arg / j, **assumptions)
    elif arg.is_constant:
        if not isinstance(expr, Superposition):
            result = expr.time(**assumptions)
        else:
            result = expr
    else:
        raise ValueError('Can only return t, f, s, omega, or jw domains')

    return result.subs(arg, **assumptions)


def transform(expr, arg, **assumptions):
    """If arg is f, s, t, omega, or jw perform domain transformation,
    otherwise perform substitution.

    Note (1 / s)(omega) will fail since 1 / s is assumed not to be
    causal and so the Fourier transform is unknown.  However,
    impedance(1 / s)(omega) will work since an impedance is assumed to
    be causal.  Alternatively, use (1 / s)(omega, causal=True). """

    arg = expr1(arg)

    new = transform1(expr, arg, **assumptions)
    return new


def call(expr, arg, **assumptions):

    if id(arg) in (id(f), id(s), id(t), id(omega), id(jw), id(jw0)):
        return expr.transform(arg, **assumptions)

    if arg in (f, s, t, omega, jw, jw0):
        return expr.transform(arg, **assumptions)

    try:
        # arg might be an int, float, complex, etc.
        if arg.has(j):
            return expr.transform(arg, **assumptions)
    except:
        pass
        
    return expr.subs(arg)


def select(expr, kind):

    if not isinstance(kind, str):
        return expr.subs(j * kind)                

    # If kind is an expr, then will add 'dc', 'time', etc. as symbols.
    
    if kind == 't':
        return expr.time()
    elif kind in ('dc', 'time'):
        return expr.subs(0)
    elif kind in ('s', 'ivp', 'super', 'laplace'):
        return expr.laplace()
    elif kind == 'f':
        return expr.fourier()
    elif kind == 'omega':
        return expr.angular_fourier()
    elif isinstance(kind, str) and kind.startswith('n'):
        return expr.angular_fourier()
    else:
        raise RuntimeError('unknown kind')


from .fexpr import FourierDomainExpression    
from .sexpr import LaplaceDomainExpression
from .texpr import TimeDomainExpression
from .omegaexpr import AngularFourierDomainExpression
from .super import Superposition
from .current import current
from .voltage import voltage
from .admittance import admittance
from .impedance import impedance
from .transfer import transfer

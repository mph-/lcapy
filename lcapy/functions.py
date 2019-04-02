import sympy as sym

def _funcwrap(func, *args):

    cls = args[0].__class__

    tweak_args = list(args)
    for m, arg in enumerate(args):
        if isinstance(arg, Expr):
            tweak_args[m] = arg.expr

    result = func(*tweak_args)

    if isinstance(args[0], Expr):
        result = cls(result)

    return result


def sin(expr):

    return _funcwrap(sym.sin, expr)


def cos(expr):

    return _funcwrap(sym.cos, expr)


def tan(expr):

    return _funcwrap(sym.tan, expr)


def atan(expr):

    return _funcwrap(sym.atan, expr)


def atan2(expr1, expr2):

    return _funcwrap(sym.atan2, expr1, expr2)


def gcd(expr1, expr2):

    return _funcwrap(sym.gcd, expr1, expr2)


def exp(expr):

    return _funcwrap(sym.exp, expr)


def sqrt(expr):

    return _funcwrap(sym.sqrt, expr)


def log(expr):

    return _funcwrap(sym.log, expr)


def log10(expr):

    return _funcwrap(sym.log, expr, 10)


def Heaviside(expr):
    """Heaviside's unit step."""

    return _funcwrap(sym.Heaviside, expr)


def H(expr):
    """Heaviside's unit step."""

    return Heaviside(expr)


def u(expr):
    """Heaviside's unit step."""

    return Heaviside(expr)


def DiracDelta(*args):
    """Dirac delta (impulse)."""

    return _funcwrap(sym.DiracDelta, *args)


def delta(expr, *args):
    """Dirac delta (impulse)."""

    return DiracDelta(expr, *args)


def conjugate(expr):
    """Complex conjugate."""
    
    return _funcwrap(sym.conjugate, expr)


from .expr import Expr

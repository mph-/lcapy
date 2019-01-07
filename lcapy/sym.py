from __future__ import division
import sympy as sym
from .sympify import canonical_name, sympify1, symbols_find
from .context import context

def sympify(expr, evaluate=True, **assumptions):
    """Create a sympy expression.

    By default, symbols are assumed to be positive unless real is
    defined.

    """
    
    if 'real' not in assumptions:
        assumptions['positive'] = True
    return sympify1(expr, context.symbols, evaluate, **assumptions)


def symbol(name, **assumptions):
    """Create a sympy symbol.

    By default, symbols are assumed to be positive unless real is
    defined.

    """
    return sympify(name, **assumptions)


def symsimplify(expr):
    """Simplify a sympy expression.  This is a hack to work around
    problems with SymPy's simplify API."""

    # Handle Matrix types
    if hasattr(expr, 'applyfunc'):
        return expr.applyfunc(lambda x: symsimplify(x))
    
    try:
        if expr.is_Function and expr.func in (sym.Heaviside, sym.DiracDelta):
            return expr
    except:
        pass

    expr = sym.simplify(expr)
    return expr

ssym = symbol('s', real=False)
tsym = symbol('t', real=True)
fsym = symbol('f', real=True)
omegasym = symbol('omega', real=True)

pi = sym.pi
j = sym.I
oo = sym.oo
inf = sym.oo

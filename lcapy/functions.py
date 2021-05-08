"""This module wraps SymPy functions and provides a few others.

Copyright 2014--2021 Michael Hayes, UCECE
"""

from .units import u as uu
import sympy as sym

class Function(object):

    def __init__(self, arg):
        self.expr = arg
    
    def __call__(self, *args):

        e_args = [expr(arg) for arg in args]
        
        cls = e_args[0].__class__

        # Handle cases like atan2(1, omega)
        if len(e_args) > 1:
            if e_args[0].is_constant:
                cls = e_args[1].__class__                
        
        result = self.expr(*[arg.expr for arg in e_args])

        result = cls(result)        
            
        if e_args[0].is_phase and self.expr in (sym.sin, sym.cos, sym.tan, sym.exp,
                                                sym.sinh, sym.cosh, sym.tanh):
            result.part = ''
        elif self.expr in (sym.atan, sym.atan2):
            result.units = uu.rad
        elif self.expr == sym.diff and isinstance(e_args[1], Expr):
            result.units = e_args[0].units / e_args[1].units
        elif self.expr == sym.integrate:
            if isinstance(e_args[1], Expr):
                result.units = e_args[0].units * e_args[1].units
            elif isinstance(e_args[1], tuple) and isinstance(e_args[1][0], Expr):
                result.units = e_args[0].units * e_args[1][0].units
        elif self.expr is sym.DiracDelta and isinstance(e_args[0], Expr):
                result.units = 1 / e_args[0].units                

        return result

    def pdb(self):
        import pdb; pdb.set_trace()
        return self
    
    
class Log10(Function):

    # TODO, figure out how to print as log10(x) rather than
    # the expansion log(x) / log(10).  This will require
    # deferment of the expansion.
    
    def __call__(self, arg):
        return super(Log10, self).__call__(arg, 10)


class SincnFunction(Function):
    """Normalised sinc function :math:`\sin(\pi x)/(\pi x)`."""

    
class SincuFunction(Function):
    """Unnormalised sinc function :math:`\sin(x)/(x)`."""

    
class PsincFunction(Function):
    """Periodic sinc function :math:`\sin(M * x)/(M * sin(x))`."""        

    
conjugate = Function(sym.conjugate)

sqrt = Function(sym.sqrt)

exp = Function(sym.exp)

log = Function(sym.log)

log10 = Log10(sym.log)

sin = Function(sym.sin)

cos = Function(sym.cos)

tan = Function(sym.tan)

cot = Function(sym.cot)

asin = Function(sym.asin)

acos = Function(sym.acos)

atan = Function(sym.atan)

atan2 = Function(sym.atan2)

acot = Function(sym.acot)

sinh = Function(sym.sinh)

cosh = Function(sym.cosh)

tanh = Function(sym.tanh)

asinh = Function(sym.asinh)

acosh = Function(sym.acosh)

atanh = Function(sym.atanh)

gcd = Function(sym.gcd)

sign = Function(sym.sign)

diff = Function(sym.diff)

integrate = Function(sym.integrate)

u = H = heaviside = Heaviside = Function(sym.Heaviside)

delta = DiracDelta = Function(sym.DiracDelta)


# def delta(expr1):
#    
#     if not hasattr(expr1, 'expr'):
#         expr1 = expr(expr1)
#     if expr1.is_discrete_time_domain or expr1.is_discrete_fourier_domain:
#         return UnitImpulse(expr1)        
#    
#     return DiracDelta(expr1)


def _ex(expr):
    if hasattr(expr, 'expr'):
        return expr.expr
    return expr


class Eq(sym.Eq):
    def __new__(cls, lhs, rhs=None, **options):
        return expr(super(Eq, cls).__new__(cls, _ex(lhs), _ex(rhs), **options))


class Add(sym.Add):
    def __new__(cls, op1, op2, **options):
        return expr(super(Add, cls).__new__(cls, _ex(op1), _ex(op2), **options))

    
class Mul(sym.Mul):
    def __new__(cls, op1, op2, **options):
        return expr(super(Mul, cls).__new__(cls, _ex(op1), _ex(op2), **options))

    
class MatAdd(sym.MatAdd):
    def __new__(cls, op1, op2, **options):
        return expr(super(MatAdd, cls).__new__(cls, _ex(op1), _ex(op2), **options))

    
class MatMul(sym.MatMul):
    def __new__(cls, op1, op2, **options):
        return expr(super(MatMul, cls).__new__(cls, _ex(op1), _ex(op2), **options))


from .extrafunctions import UnitImpulse, UnitStep
from .extrafunctions import sincn as sincn1
from .extrafunctions import sincu as sincu1
from .extrafunctions import psinc as psinc1
from .extrafunctions import rect as rect1
from .extrafunctions import tri as tri1
from .extrafunctions import trap as trap1


ui = unitimpulse = Function(UnitImpulse)

us = unitstep = Function(UnitStep)

rect = Function(rect1)

sinc = SincnFunction(sincn1)

sincn = SincnFunction(sincn1)

sincu = SincuFunction(sincu1)

psinc = PsincFunction(psinc1)

tri = Function(tri1)

trap = Function(trap1)
 

from .expr import Expr, expr

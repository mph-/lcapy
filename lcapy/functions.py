"""This module wraps SymPy functions and provides a few others.
The underlying SymPy function is obtained with the sympy attribute.

Copyright 2014--2022 Michael Hayes, UCECE
"""

from .units import u as uu
import sympy as sym

__all__ = ('conjugate', 'sqrt', 'exp', 'log', 'log10', 'sin', 'cos',
           'tan', 'cot', 'asin', 'acos', 'atan', 'atan2', 'acot',
           'sinh', 'cosh', 'tanh', 'coth', 'asinh', 'acosh', 'atanh',
           'acoth', 'gcd', 'sign', 'diff', 'integrate', 'u', 'H',
           'heaviside', 'Heaviside', 'delta', 'DiracDelta', 'ui',
           'us', 'unitimpulse', 'unitstep', 'UnitImpulse', 'UnitStep',
           'rect', 'sinc', 'sincn', 'sincu', 'psinc', 'tri', 'trap',
           'ramp', 'rampstep', 'Sum', 'dtrect', 'dtsign', 'Piecewise',
           'Eq', 'Ne', 'Lt', 'Le', 'Gt', 'Ge', 'Derivative',
           'Integral', 'Max', 'Min', 're', 'im', 'MatMul', 'MatAdd',
           'degrees', 'radians', 'limit')

from .extrafunctions import Degrees as Degrees1
from .extrafunctions import Radians as Radians1
from .extrafunctions import UnitImpulse as UnitImpulse1
from .extrafunctions import UnitStep as UnitStep1
from .extrafunctions import sincn as sincn1
from .extrafunctions import sincu as sincu1
from .extrafunctions import psinc as psinc1
from .extrafunctions import rect as rect1
from .extrafunctions import dtrect as dtrect1
from .extrafunctions import dtsign as dtsign1
from .extrafunctions import tri as tri1
from .extrafunctions import trap as trap1
from .extrafunctions import ramp as ramp1
from .extrafunctions import rampstep as rampstep1


class Function(object):

    def __init__(self, arg):

        if isinstance(arg, str):
            # Handled AppliedUndef
            arg = sym.Function(arg)
        self.expr = arg

    def __repr__(self):

        return repr(self.expr)

    @property
    def sympy(self):
        """Return SymPy function."""

        return self.expr

    def __call__(self, *args, **kwargs):

        func = self.expr

        e_args = list(map(expr, args))

        if func == sym.Piecewise:
            cls = e_args[0][0].__class__
        elif func in (sym.Integral, sym.integrate):
            # Integral(const, (tau, 0, t)) -> use t class
            # Integral(integrand, (tau, a, b)) -> use integrand class
            if e_args[0].is_constant:
                cls = e_args[1][2].__class__
            else:
                cls = e_args[0].__class__
        else:
            cls = e_args[0].__class__

            # Handle cases like atan2(1, omega)
            if len(e_args) > 1:
                if e_args[0].is_constant:
                    cls = e_args[1].__class__

        try:
            if cls.is_discrete_time_domain or cls.is_discrete_fourier_domain:

                # TODO: fix for Heaviside change
                for old, new in function_mapping.items():
                    if func is old:
                        func = new
                        break
        except:
            pass

        s_args = delcapify(e_args)
        result = func(*s_args, **kwargs)

        result = cls(result, **e_args[0].assumptions)

        if (self.expr in (sym.sin, sym.cos, sym.tan, sym.cot, sym.exp,
                          sym.sinh, sym.cosh, sym.tanh, sym.coth)
                and e_args[0].is_phase):
            result.part = ''
        elif self.expr in (sym.atan, sym.atan2):
            result.units = uu.rad
        elif self.expr == sym.diff and isinstance(e_args[1], Expr):
            result.units = e_args[0].units / e_args[1].units
        elif self.expr == Degrees1:
            result.units = uu.deg
        elif self.expr == Radians1:
            result.units = uu.rad
        elif self.expr == sym.integrate:
            if isinstance(e_args[1], Expr):
                result.units = e_args[0].units * e_args[1].units
            elif isinstance(e_args[1], tuple) and isinstance(e_args[1][0], Expr):
                result.units = e_args[0].units * e_args[1][0].units
        elif self.expr is sym.DiracDelta and isinstance(e_args[0], Expr):
            result.units = 1 / e_args[0].units

        return result

    def pdb(self):
        import pdb
        pdb.set_trace()
        return self


class Log10(Function):

    # TODO, figure out how to print as log10(x) rather than
    # the expansion log(x) / log(10).  This will require
    # deferment of the expansion.

    def __call__(self, arg):
        return super(Log10, self).__call__(arg, 10)


def function_wrap(symfunc):
    """Create an Lcapy function from a SymPy function.

    An Lcapy function works with Lcapy expressions and returns an
    Lcapy expression.  It acquires the SymPy docstring.

    """

    name = symfunc.__name__
    docstring = symfunc.__doc__.replace('sympy', 'lcapy')

    newclass = type(name, (Function, ), {'__doc__': docstring})
    return newclass(symfunc)


conjugate = function_wrap(sym.conjugate)

sqrt = function_wrap(sym.sqrt)

exp = function_wrap(sym.exp)

log = function_wrap(sym.log)

log10 = Log10(sym.log)

sin = function_wrap(sym.sin)

cos = function_wrap(sym.cos)

tan = function_wrap(sym.tan)

cot = function_wrap(sym.cot)

asin = function_wrap(sym.asin)

acos = function_wrap(sym.acos)

atan = function_wrap(sym.atan)

atan2 = function_wrap(sym.atan2)

acot = function_wrap(sym.acot)

sinh = function_wrap(sym.sinh)

cosh = function_wrap(sym.cosh)

tanh = function_wrap(sym.tanh)

coth = function_wrap(sym.coth)

asinh = function_wrap(sym.asinh)

acosh = function_wrap(sym.acosh)

atanh = function_wrap(sym.atanh)

acoth = function_wrap(sym.acoth)

gcd = function_wrap(sym.gcd)

sign = function_wrap(sym.sign)

degrees = function_wrap(Degrees1)

radians = function_wrap(Radians1)

diff = function_wrap(sym.diff)

integrate = function_wrap(sym.integrate)

limit = function_wrap(sym.limit)

u = H = heaviside = Heaviside = function_wrap(sym.Heaviside)

delta = DiracDelta = function_wrap(sym.DiracDelta)

Piecewise = function_wrap(sym.Piecewise)

Derivative = function_wrap(sym.Derivative)

Integral = function_wrap(sym.Integral)

Eq = function_wrap(sym.Eq)

Ne = function_wrap(sym.Ne)

Lt = function_wrap(sym.Lt)

Le = function_wrap(sym.Le)

Gt = function_wrap(sym.Gt)

Ge = function_wrap(sym.Ge)

Add = function_wrap(sym.Add)

Mul = function_wrap(sym.Mul)

MatAdd = function_wrap(sym.MatAdd)

MatMul = function_wrap(sym.MatMul)

Sum = function_wrap(sym.Sum)

Min = function_wrap(sym.Min)

Max = function_wrap(sym.Max)

re = function_wrap(sym.re)

im = function_wrap(sym.im)

MatMul = function_wrap(sym.MatMul)

MatAdd = function_wrap(sym.MatAdd)

ui = unitimpulse = UnitImpulse = function_wrap(UnitImpulse1)

us = unitstep = UnitStep = function_wrap(UnitStep1)

rect = function_wrap(rect1)

dtrect = function_wrap(dtrect1)

dtsign = function_wrap(dtsign1)

sinc = sincn = function_wrap(sincn1)

sincu = function_wrap(sincu1)

psinc = function_wrap(psinc1)

tri = function_wrap(tri1)

trap = function_wrap(trap1)

ramp = function_wrap(ramp1)

rampstep = function_wrap(rampstep1)

function_mapping = {sym.Heaviside: UnitStep1,
                    sym.DiracDelta: UnitImpulse1,
                    rect1: dtrect,
                    sym.sign: dtsign}

from .expr import Expr, expr, delcapify  # nopep8

"""
Copyright 2014, 2015, 2016 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from lcapy.core import t, s, Vs, Is, Zs, Ys, cExpr, sExpr, tExpr, tsExpr, cos, exp, symbol, j, Vphasor, Iphasor, Yphasor, Zphasor, omega1

class Network(object):

    voltage_source = False
    current_source = False

    # True if initial conditions are zero (or unspecified).
    zeroic = True

    # None if component does not have initial conditions.
    # True if initial conditions are specified.
    # False if initial conditions are not specified.
    hasic = None

    def __init__(self, args):

        if not hasattr(self, 'args'):
            self.args = args

    def _tweak_args(self):

        if not hasattr(self, 'args'):
            return ()

        args = self.args
        # Drop the initial condition for L or C if it is zero.
        if isinstance(self, (L, C)) and args[1] == 0:
            args = args[:-1]

        modargs = []
        for arg in args:
            if isinstance(arg, sExpr):
                arg = arg.expr

            modargs.append(arg)
        return modargs

    def __repr__(self):

        argsrepr = ', '.join([arg.__repr__() for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def __str__(self):

        def fmt(arg):
            if False and isinstance(arg, str):
                return "'" + arg + "'"
            return arg.__str__()

        argsrepr = ', '.join([fmt(arg) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def _repr_pretty_(self, p, cycle):

        p.text(self.pretty())

    def _repr_latex_(self):

        return '$%s$' % self.latex()

    def pretty(self):

        argsrepr = ', '.join([sym.pretty(arg) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def latex(self):

        argsrepr = ', '.join([latex_str(sym.latex(arg)) for arg in self._tweak_args()])
        return '\\mathrm{%s}(%s)' % (self.__class__.__name__, argsrepr)

    def simplify(self):

        return self

    @property
    def Vac(self):
        return 0

    @property
    def Iac(self):
        return 0

    @property
    def Yac(self):
        return self.Y.jomega()

    @property
    def Zac(self):
        return self.Z.jomega()

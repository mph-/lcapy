"""
Copyright 2014, 2015, 2016 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from lcapy.core import t, s, Vs, Is, Zs, Ys, cExpr, sExpr, tExpr, tsExpr, cos, exp, symbol, j, Vphasor, Iphasor, Yphasor, Zphasor, omega1
from lcapy.schematic import Schematic

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
        self.node_counter = 0


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

    @property
    def elements(self):
        raise ValueError('TODO')

    @property 
    def node(self):

        self.node_counter += 1
        return self.node_counter

    def netargs(self):

        def quote(arg):

            if ('(' in arg) or (')' in arg) or (' ' in arg) or (',' in arg):
                return '{%s}' % arg
            return arg

        return ' '.join([quote(str(arg)) for arg in self.args])

    def netlist(self, n1=None, n2=None):

        if n1 == None:
            n1 = self.node
        if n2 == None:
            n2 = self.node

        netname = self.__class__.__name__ if self.netname == '' else self.netname
        return '%s %s %s %s %s; right' % (netname, n1, n2, 
                                          self.netkeyword, self.netargs())

    @property
    def sch(self):

        if hasattr(self, '_sch'):
            return self._sch

        netlist = self.netlist()
        sch = Schematic()
        for net in netlist.split('\n'):
            sch.add(net)
        self._sch = sch
        return sch

    def draw(self, filename=None, label_ids=False,
             label_values=True, draw_nodes='connections',
             label_nodes=False):

        self.node_counter = 0
        self.sch.draw(filename=filename, label_ids=label_ids, 
                      label_values=label_values, 
                      draw_nodes=draw_nodes, label_nodes=label_nodes)
        
    @property
    def cct(self):

        if hasattr(self, '_cct'):
            return self._cct

        from lcapy.netlist import Circuit

        self.node_counter = 0
        netlist = self.netlist(self)
        cct = Circuit()
        for net in netlist.split('\n'):
            cct.add(net)

        # Create ground reference.
        cct.add('W %d 0' % (self.node - 1))
        self._cct = cct
        return cct


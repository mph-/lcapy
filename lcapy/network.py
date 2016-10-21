"""
Copyright 2014, 2015, 2016 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from lcapy.core import Vphasor, Iphasor, sExpr
from lcapy.latex import latex_str
from lcapy.schematic import Schematic
from lcapy.netlist import Netlist
from lcapy.circuit import Circuit

class Network(Netlist):

    voltage_source = False
    current_source = False

    # True if initial conditions are zero (or unspecified).
    zeroic = True

    # None if component does not have initial conditions.
    # True if initial conditions are specified.
    # False if initial conditions are not specified.
    hasic = None

    def __init__(self):

        super(Network, self).__init__()

    def _tweak_args(self):

        if not hasattr(self, 'args'):
            return ()

        modargs = []
        for arg in self.args:
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
    def Vocac(self):
        return Vphasor(0)

    @property
    def Iscac(self):
        return Iphasor(0)

    @property
    def Yac(self):
        return self.Y.jomega()

    @property
    def Zac(self):
        return self.Z.jomega()

    def _add_elements(self):

        netlist = self.netlist()
        for net in netlist.split('\n'):
            self._add(net)

        # Hack, create ground reference.
        self._add('W %d 0' % (self.node - 1))

    @property 
    def node(self):

        if not hasattr(self, 'node_counter'):
            self.node_counter = 0
        
        self.node_counter += 1
        return self.node_counter

    def netargs(self):

        def quote(arg):

            if ('(' in arg) or (')' in arg) or (' ' in arg) or (',' in arg):
                return '{%s}' % arg
            return arg

        return ' '.join([quote(str(arg)) for arg in self.args])

    def net_make(self, net, n1=None, n2=None):

        if n1 == None:
            n1 = net.node
        if n2 == None:
            n2 = net.node

        netname = self.__class__.__name__ if self.netname == '' else self.netname
        return '%s %s %s %s %s; right' % (netname, n1, n2, 
                                          self.netkeyword, self.netargs())

    def netlist(self):

        self.node_counter = 0
        return self.net_make(self)


    @property
    def sch(self):
        """Convert a Network object into a Schematic object."""

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

        self.sch.draw(filename=filename, label_ids=label_ids, 
                      label_values=label_values, 
                      draw_nodes=draw_nodes, label_nodes=label_nodes)
        
    @property
    def cct(self):
        """Convert a Network object into a Circuit object for debugging."""

        if hasattr(self, '_cct'):
            return self._cct

        netlist = self.netlist()
        cct = Circuit()
        for net in netlist.split('\n'):
            cct.add(net)
        self._cct = cct
        return cct

"""This module provides support for the common aspects of Circuit and
Network classes.

Copyright 2014, 2015 Michael Hayes, UCECE

"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from lcapy.core import pprint, Hs, Vs, Zs, Ys, Expr, tsym
from lcapy.core import s, j, omega, uppercase_name, global_context
from lcapy.schematic import Schematic, Opts, SchematicOpts
from lcapy.mna import MNA
from lcapy.netfile import NetfileMixin
import lcapy.mnacpts as cpts
import re
from copy import copy

class Ldict(dict):

    """Lazy dictionary for inverse Laplace."""

    def __init__(self, Vdict, **assumptions):

        super(Ldict, self).__init__()

        self.Vdict = Vdict
        self.assumptions = assumptions

    def __getitem__(self, key):

        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key
        
        # Note, need to use keys method to catch branch names.
        if (key not in self) and (key in self.Vdict.keys()):
            v = self.Vdict[key].time(**self.assumptions)
            self[key] = v
            return v

        return super(Ldict, self).__getitem__(key)

    def keys(self):

        return list(self.Vdict)


class Node(object):

    def __init__(self, cct, name):

        self.cct = cct
        self.name = name
        self.pos = None
        self.port = False
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        # List of elements connected to this node.
        self.list = []

    @property
    def V(self):
        """Node voltage with respect to ground."""

        self.cct._solve()
        return self.cct._V[self.name]

    @property
    def v(self):
        """Node time-domain voltage with respect to ground."""

        return self.cct.v[self.name]

    def append(self, cpt):

        if cpt.type in ('P', ):
            self.port = True

        self.list.append(cpt)


class Netlist(MNA, NetfileMixin):

    def __init__(self, filename=None):

        self._elements = {}
        self.nodes = {}
        # Shared nodes (with same voltage)
        self.snodes = {}
        self.context = global_context.new()
        self._init_parser(cpts)

        self.opts = SchematicOpts()

        if filename is not None:
            self.netfile_add(filename)

    def __getitem__(self, name):
        """Return element or node by name."""

        # If name is an integer, convert to a string.
        if isinstance(name, int):
            name = '%d' % name

        if name in self.nodes:
            return self.nodes[name]

        if name in self._elements:
            return self._elements[name]

        # Try first anonymous name.
        if name + 'anon1' in self._elements:
            return self._elements[name + 'anon1']

        raise AttributeError('Unknown element or node name %s' % name)

    def __getattr__(self, attr):
        """Return element or node by name."""

        # This gets called if there is no explicit attribute attr for
        # this instance.  This is primarily for accessing elements
        # and non-numerical node names.

        return self.__getitem__(attr)

    def __repr__(self):
        
        return self.netlist()

    @property
    def elements(self):

        if hasattr(self, '_add_elements'):
            if self._elements == {}:
                self._add_elements()

        return self._elements

    def netlist(self):
        """Return the current netlist."""

        return '\n'.join([str(cpt) for cpt in self._elements.values()])

    def _node_add(self, node, cpt):

        if node not in self.nodes:
            self.nodes[node] = Node(self, node)
        self.nodes[node].append(cpt)

        vnode = self.nodes[node].rootname

        if vnode not in self.snodes:
            self.snodes[vnode] = []

        if node not in self.snodes[vnode]:
            self.snodes[vnode].append(node)

    def _invalidate(self):

        for attr in ('_sch', '_A', '_V', '_I', '_Vd', '_node_list',
                     '_vcache', '_icache'):
            if hasattr(self, attr):
                delattr(self, attr)

    def _cpt_add(self, cpt):

        self._invalidate()

        opts = Opts(cpt.opts_string)
        cpt.opts = opts

        if cpt.name in self._elements:
            print('Overriding component %s' % cpt.name)
            # Need to search lists and update component.
            # For example, remove nodes that are only connected
            # to this component.
        else:
            # Check that this name won't conflict with an attr.
            # For example, cannot have name V or I.  Perhaps
            # rename these attributes?
            if hasattr(self, cpt.name):
                raise ValueError('Invalid component name %s' % cpt.name)

        self._elements[cpt.name] = cpt

        for node in cpt.nodes:
            self._node_add(node, cpt)

    def remove(self, name):
        """Remove specified element."""

        self._invalidate()

        if name not in self._elements:
            raise ValueError('Unknown component: ' + name)
        self._elements.pop(name)
        # TODO, remove nodes that are only connected
        # to this component.

    @property
    def v(self):
        """Return dictionary of t-domain node voltages indexed by node name
        and voltage differences indexed by branch name."""

        if not hasattr(self, '_vcache'):
            self._vcache = Ldict(self.V, **self.assumptions)

        return self._vcache

    @property
    def i(self):
        """Return dictionary of t-domain branch currents indexed
        by component name."""

        if not hasattr(self, '_icache'):
            self._icache = Ldict(self.I, **self.assumptions)

        return self._icache

    def Voc(self, Np, Nm):
        """Return open-circuit s-domain voltage between nodes Np and Nm."""

        return self.V[Np] - self.V[Nm]

    def voc(self, Np, Nm):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).time()

    def Isc(self, Np, Nm):
        """Return short-circuit s-domain current between nodes Np and Nm."""

        self._add('Vshort_ %d %d 0' % (Np, Nm))

        Isc = self.Vshort_.I
        self.remove('Vshort_')

        return Isc

    def isc(self, Np, Nm):
        """Return short-circuit t-domain current between nodes Np and Nm."""

        return self.Isc(Np, Nm).time()

    def thevenin(self, Np, Nm):
        """Return Thevenin model between nodes Np and Nm."""

        from lcapy.oneport import V, Z

        Voc = self.Voc(Np, Nm)

        return V(Voc) + Z(self.impedance(Np, Nm))

    def norton(self, Np, Nm):
        """Return Norton model between nodes Np and Nm."""

        from lcapy.oneport import I, Y

        Isc = self.Isc(Np, Nm)

        return I(Isc) | Y(self.admittance(Np, Nm))

    def admittance(self, Np, Nm):
        """Return admittance between nodes Np and Nm with independent 
        sources killed.

        """

        new = self.kill()

        # Connect 1 V s-domain voltage source between nodes and
        # measure current.
        new._add('Vin_ %d %d {2 * DiracDelta(t)}' % (Np, Nm))
        If = -new.Vin_.I
        new.remove('Vin_')

        Y = Ys(If, causal=True)
        return Y

    def impedance(self, Np, Nm):
        """Return impedance between nodes Np and Nm with independent
        sources killed.

        """

        new = self.kill()

        # Connect 1 A s-domain current source between nodes and
        # measure voltage.
        new._add('Iin_ %d %d {2 * DiracDelta(t)}' % (Np, Nm))
        Vf = new.Voc(Np, Nm)
        new.remove('Iin_')

        Z = Zs(Vf, causal=True)
        return Z

    def transfer(self, N1p, N1m, N2p, N2m):
        """Create voltage transfer function V2 / V1 where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed."""

        new = self.kill()
        new._add('V1_ %d %d s 1' % (N1p, N1m))

        H = Hs(new.Voc(N2p, N2m) / new.V1_.V, causal=True)

        return H

    def Amatrix(self, N1p, N1m, N2p, N2m):
        """Create A matrix from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """

        from lcapy.twoport import AMatrix

        if self.Voc(N1p, N1m) != 0 or self.Voc(N2p, N2m) != 0:
            raise ValueError('Network contains independent sources')

        try:

            self._add('V1_ %d %d s 1' % (N1p, N1m))

            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = Hs(self.V1_.V / self.Voc(N2p, N2m))

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = Zs(self.V1_.V / self.Isc(N2p, N2m))

            self.remove('V1_')

            self._add('I1_ %d %d s 1' % (N1p, N1m))

            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure I2 with port 2 open-circuit
            A21 = Ys(-self.I['I1_'] / self.Voc(N2p, N2m))

            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure I2 with port 2 short-circuit
            A22 = Hs(-self.I['I1_'] / self.Isc(N2p, N2m))

            self.remove('I1_')
            return AMatrix(A11, A12, A21, A22)

        except ValueError:
            raise ValueError('Cannot create A matrix')

    def _kill(self, sourcenames):

        new = Netlist()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            if cpt.name in sourcenames:
                net = cpt.kill()
            else:
                net = cpt.kill_initial()
            new._add(net)
        return new        

    def kill_except(self, *args):
        """Return a new circuit with all but the specified sources killed;
        i.e., make the voltage sources short-circuits and the current
        sources open-circuits.  If no sources are specified, all
        independent sources (including initial conditions) are killed.

        """

        for arg in args:
            if arg not in self.independent_sources:
                raise ValueError('Element %s is not a known source' % arg)
        sources = []
        for cpt in self.independent_sources.values():
            if cpt.name not in args:
                sources.append(cpt.name)
        return self._kill(sources)

    def kill(self, *args):
        """Return a new circuit with the specified sources killed; i.e., make
        the voltage sources short-circuits and the current sources
        open-circuits.  If no sources are specified, all independent
        sources (including initial conditions) are killed.

        """

        if len(args) == 0:
            return self.kill_except()

        sources = []
        for arg in args:
            if arg not in self.independent_sources:
                raise ValueError('Element %s is not a known source' % arg)
            sources.append(self.independent_sources[arg].name)

        return self._kill(sources)

    def twoport(self, N1p, N1m, N2p, N2m):
        """Create twoport model from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """

        from lcapy.twoport import TwoPortBModel

        V2b = self.Voc(N2p, N2m)
        I2b = self.Isc(N2p, N2m)

        A = self.kill().Amatrix(N1p, N1m, N2p, N2m)

        return TwoPortBModel(A.B, V2b, I2b)

    @property
    def sch(self):

        if hasattr(self, '_sch'):
            return self._sch

        sch = Schematic()

        netlist = self.netlist()
        for net in netlist.split('\n'):
            sch.add(net)

        self._sch = sch
        return sch

    def pre_initial_model(self):
        """Generate circuit model for determining the pre-initial
        conditions."""

        new = Netlist()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            net = cpt.pre_initial_model()
            new._add(net)
        return new        

    def s_model(self, var=s):

        new = Netlist()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            net = cpt.s_model(var)
            new._add(net)
        return new        

    def ac_model(self):
        return self.s_model(j * omega)

    def draw(self, filename=None, **kwargs):

        cct = self
        if kwargs.pop('s_model', False):
            cct = cct.s_model()

        return cct.sch.draw(filename=filename, opts=self.opts, **kwargs)

    @property
    def is_causal(self):
        """Return True if all independent sources are causal and not an
        initial value problem (unless all the initial values are zero)."""

        # If some of the initial conditions are specified and are non-zero
        # then have a non-causal initial value problem.
        if not self.zeroic:
            return False

        independent_sources = self.independent_sources
        for cpt in self.independent_sources.values():
            if not cpt.is_causal:
                return False
        return True

    @property
    def is_dc(self):
        """Return True if all independent sources are DC and not an
        initial value problem.  The initial value problem may collapse
        to a DC problem but we cannot prove this yet."""

        if self.initial_value_problem:
            return False

        independent_sources = self.independent_sources
        for cpt in independent_sources.values():
            if not cpt.is_dc:
                return False
        return True

    @property
    def is_ac(self):
        """Return True if all independent sources are AC and not an
        initial value problem."""

        if self.initial_value_problem:
            return False

        independent_sources = self.independent_sources
        for cpt in independent_sources.values():
            if not cpt.is_ac:
                return False
        return True

    @property
    def assumptions(self):

        assumptions = {}
        if self.is_ac:
            assumptions['ac'] = True
        if self.is_dc:
            assumptions['dc'] = True
        if self.is_causal:
            assumptions['causal'] = True
        return assumptions

    @property
    def zeroic(self):
        """Return True if the initial conditions for all components are zero."""

        for cpt in self.elements.values():
            if not cpt.zeroic:
                return False
        return True

    @property
    def initial_value_problem(self):
        """Return True if any components that allow initial conditions
        have them explicitly defined."""

        for cpt in self.elements.values():
            if cpt.hasic is None:
                continue
            if cpt.hasic:
                return True

        return False

    @property
    def missing_ic(self):
        """Return components that allow initial conditions but do not have
        them explicitly defined."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.hasic is False)

    @property
    def explicit_ic(self):
        """Return components that have explicitly defined initial conditions."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.hasic is True)

    @property
    def allow_ic(self):
        """Return components (L and C) that allow initial conditions."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.hasic is not None)

    @property
    def noncausal_sources(self):
        """Return non-causal independent sources."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.source and not cpt.is_causal)

    @property
    def independent_sources(self):
        """Return independent sources (this does not include
        implicit sources due to initial conditions)."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.source)


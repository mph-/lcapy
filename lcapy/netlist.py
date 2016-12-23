"""This module provides support for the common aspects of Circuit and
Network classes.

Copyright 2014-2016 Michael Hayes, UCECE

"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from lcapy.core import pprint, Hs, Vs, Zs, Ys, Expr, tsym, Vt, It
from lcapy.core import s, j, omega, uppercase_name, global_context
from lcapy.core import Vtype, sqrt
from lcapy.schematic import Schematic, Opts, SchematicOpts
from lcapy.mna import MNA
from lcapy.netfile import NetfileMixin
import lcapy.mnacpts as cpts
import re
from copy import copy
from collections import OrderedDict


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

        return self.cct.get_Vd(self.name, '0')

    @property
    def v(self):
        """Node time-domain voltage with respect to ground."""

        return self.cct.get_vd(self.name, '0')

    def append(self, cpt):

        if cpt.type in ('P', ):
            self.port = True

        self.list.append(cpt)


class NetlistMixin(object):

    def __init__(self, filename=None, context=None):

        self._elements = OrderedDict()
        self.nodes = {}
        if context is None:
            context = global_context.new()
        
        self.context = context
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
        """Return element or node by name.  This gets called if there is no
        explicit attribute attr for this instance.  This is primarily
        for accessing elements and non-numerical node names.  It also
        gets called if the called attr throws an AttributeError
        exception.  The annoying thing is that hasattr uses getattr
        and checks for an exception."""

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

    def _invalidate(self):

        for attr in ('_sch', '_A', '_V', '_I', '_node_list'):
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

    def copy(self):
        """Create a copy of the netlist"""

        new = self._new()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            new._add(str(cpt))
        return new        

    def _new(self):

        # TODO.  Copy or share?
        context = self.context
        if self.__class__ == 'Circuit':
            return Circuit(context=context)
        # If have OnePort, Network, etc., treat as Netlist
        return Netlist(context=context)

    def remove(self, name):
        """Remove specified element."""

        self._invalidate()

        if name not in self._elements:
            raise ValueError('Unknown component: ' + name)
        self._elements.pop(name)
        # TODO, remove nodes that are only connected
        # to this component.


    def Voc(self, Np, Nm):
        """Return open-circuit s-domain voltage between nodes Np and Nm."""

        return Vtype(self.V[Np] - self.V[Nm],  **self.assumptions)

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
        new._add('Vin_ %d %d {DiracDelta(t)}' % (Np, Nm))
        If = -new.Vin_.I
        new.remove('Vin_')

        return Ys(If)

    def impedance(self, Np, Nm):
        """Return impedance between nodes Np and Nm with independent
        sources killed.

        """

        new = self.kill()

        # Connect 1 A s-domain current source between nodes and
        # measure voltage.
        new._add('Iin_ %d %d {DiracDelta(t)}' % (Np, Nm))
        Vf = new.Voc(Np, Nm)
        new.remove('Iin_')

        return Zs(Vf)

    def transfer(self, N1p, N1m, N2p, N2m):
        """Create voltage transfer function V2 / V1 where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed."""

        new = self.kill()
        new._add('V1_ %d %d {DiracDelta(t)}' % (N1p, N1m))

        V2 = new.Voc(N2p, N2m)
        V1 = new.V1_.V

        return Hs(V2 / V1, causal=True)

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
            self._add('V1_ %d %d {DiracDelta(t)}' % (N1p, N1m))

            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = Hs(self.V1_.V / self.Voc(N2p, N2m))

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = Zs(self.V1_.V / self.Isc(N2p, N2m))

            self.remove('V1_')

            self._add('I1_ %d %d {DiracDelta(t)}' % (N1p, N1m))

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

        new = self._new()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            if cpt.name in self.control_sources:
                net = cpt.zero()                
            elif cpt.name in sourcenames:
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

    def _noisy(self, resistornames):

        new = self._new()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            if cpt.name in resistornames:
                net = cpt.noisy()
            else:
                net = str(cpt)
            new._add(net)
        return new        

    def noisy_except(self, *args):
        """Return a new circuit with all but the specified resistors in series
        with noise voltage sources"""

        for arg in args:
            if arg not in self.elements and self.elements[arg].type != 'R':
                raise ValueError('Element %s is not a known resistor' % arg)
        resistors = []
        for cpt in self.elements.values():
            if cpt.type == 'R' and cpt.name not in args:
                resistors.append(cpt.name)
        return self._noisy(resistors)

    def noisy(self, *args):
        """Return a new circuit with the specified resistors in series
        with noise voltage sources"""

        if len(args) == 0:
            return self.noisy_except()

        resistors = []
        for arg in args:
            if arg not in self.elements and self.elements[arg].type != 'R':
                raise ValueError('Element %s is not a known resistor' % arg)
            resistors.append(arg)

        return self._noisy(resistors)
    

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

        new = self._new()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            net = cpt.pre_initial_model()
            new._add(net)
        return new        

    def s_model(self, var=s):

        new = self._new()
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

        if self.is_ivp:
            return False

        independent_sources = self.independent_sources
        if independent_sources == {}:
            return False

        for cpt in independent_sources.values():
            if not cpt.is_dc:
                return False
        return True

    @property
    def is_ac(self):
        """Return True if all independent sources are AC and not an
        initial value problem."""

        if self.is_ivp:
            return False

        independent_sources = self.independent_sources
        if independent_sources == {}:
            return False

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
    def is_ivp(self):
        """Return True for an initial value problem.  This is True if any
        component that allows initial conditions has them explicitly
        defined.

        """
        return self.initial_value_problem

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
        """Return dictionary of non-causal independent sources."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.source and not cpt.is_causal)

    @property
    def independent_sources(self):
        """Return dictionary of independent sources (this does not include
        implicit sources due to initial conditions)."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.source)

    @property
    def control_sources(self):
        """Return dictionary of voltage sources required to specify control
        current for CCVS and CCCS components."""

        result = {}

        for key, cpt in self.elements.items():
            if cpt.need_control_current:
                result[cpt.args[0]] = self.elements[cpt.args[0]]
        return result

    def sources_kind(self, kind):
        """Return sources that match kind"""

        result = []
        for source, cpt in self.independent_sources.items():
            if cpt.is_dc and kind == 'dc':
                result.append(source)
            elif cpt.is_ac and kind == 'ac':
                result.append(source)
            elif cpt.is_noisy and kind == 'n':
                result.append(source)
            elif not cpt.is_dc and not cpt.is_ac and kind == 's':
                result.append(source)
        return result


class Domains(dict):

    def __getattr__(self, attr):

        if attr not in self:
            raise AttributeError('Unknown attribute %s' % attr)
        return self[attr]
    
                        
class Netlist(NetlistMixin, NetfileMixin):

    def __init__(self, filename=None, context=None):

        super (Netlist, self).__init__(filename, context)
        self._invalidate()

    def _invalidate(self):
        self.sub = Domains({'s' : None, 'dc' : None,
                            'ac' : None, 'n' : None})

    def _sub_make(self, kind):

        result = {}
        sources = self.sources_kind(kind)
        for source in sources:
            result[source] = SubNetlist(self, source, kind)
        return result

    def _subnetlist(self, kind):
        if self.sub[kind] is None:
            self.sub[kind] = self._sub_make(kind)
        return self.sub[kind]        

    @property
    def dc(self):
        return self._subnetlist('dc')

    @property        
    def s(self):
        return self._subnetlist('s')        

    @property        
    def ac(self):
        return self._subnetlist('ac')

    @property        
    def n(self):
        return self._subnetlist('n')

    @property
    def V(self):
        """Return dictionary of node voltages for each frequency domain"""

        # Hack to generate entries in sub
        self.ac
        self.dc
        self.s
        self.n

        result = {}

        for kind, sub in self.sub.items():
            for source, subnetlist in sub.items():
                for node, value in subnetlist.V.items():
                    if node not in result:
                        result[node] = Domains({'s' : 0, 'dc' : 0, 'ac' : 0, 'n' : 0})
                    if kind == 'n':
                        result[node][kind] += value**2
                    else:
                        result[node][kind] += value

        for node, domains in result.items():                        
            domains['n'] = sqrt(domains['n'])
        return result

    def get_I(self, name):
        """Current through component"""

        self.ac
        self.dc
        self.s
        self.n
        result = Domains({'s' : 0, 'dc' : 0, 'ac' : 0, 'n' : 0})        

        for kind, sub in self.sub.items():
            for source, subnetlist in sub.items():
                I = subnetlist.get_I(name)
                if kind == 'n':
                    I = I * I
                result[kind] += I

        result['n'] = sqrt(result['n'])
        return result        

    def get_i(self):
        """Time-domain current through component"""

        I = self.get_I(n2, n1)

        # TODO, integrate noise
        result = 0
        if I.s != 0:
            result += I.s.time()
        if I.ac != 0:
            result += I.ac.time()
        if I.dc != 0:
            result += I.dc.time()
        return It(result)

    def get_Vd(self, n2, n1):
        """Voltage drop between nodes"""

        self.ac
        self.dc
        self.s
        self.n        
        result = Domains({'s' : 0, 'dc' : 0, 'ac' : 0, 'n' : 0})        

        for kind, sub in self.sub.items():
            for source, subnetlist in sub.items():
                Vd = subnetlist.get_Vd(n2, n1)
                if kind == 'n':
                    Vd = Vd * Vd
                result[kind] += Vd

        result['n'] = sqrt(result['n'])
        return result        

    def get_vd(self, n2, n1):
        """Time-domain voltage drop between nodes"""

        Vd = self.get_Vd(n2, n1)

        # TODO, integrate noise
        result = 0
        if Vd.s != 0:
            result += Vd.s.time()
        if Vd.ac != 0:
            result += Vd.ac.time()
        if Vd.dc != 0:
            result += Vd.dc.time()
        return Vt(result)        

    
class SubNetlist(NetlistMixin, MNA):

    def __new__(cls, netlist, sources, kind):

        obj = netlist.kill_except(sources)
        # Hack
        obj.__class__ = cls
        return obj

    def __init__(cls, netlist, sources, kind):
        pass

    def get_I(self, name):
        """Current through component"""

        self._solve()
        return self._I[name]

    def get_i(self):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def get_Vd(self, n2, n1):
        """Voltage drop between nodes"""

        self._solve()
        return self._V[n2] - self._V[n1]

    def get_vd(self, n2, n1):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(n2, n1).time()

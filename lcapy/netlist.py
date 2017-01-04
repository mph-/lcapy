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
from lcapy.core import Vsuper, Isuper
from lcapy.schematic import Schematic, Opts, SchematicOpts
from lcapy.mna import MNA, Nodedict, Branchdict
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

    def _cpt_add(self, cpt):

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

        return self.get_Vd(Np, Nm)

    def voc(self, Np, Nm):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).time()

    def Isc(self, Np, Nm):
        """Return short-circuit s-domain current between nodes Np and Nm."""

        self.add('Vshort_ %d %d 0' % (Np, Nm))

        Isc = self.Vshort_.I * Ys(1)
        self.remove('Vshort_')

        return Isc

    def isc(self, Np, Nm):
        """Return short-circuit t-domain current between nodes Np and Nm."""

        return self.Isc(Np, Nm).time()

    def thevenin(self, Np, Nm):
        """Return s-domain Thevenin model between nodes Np and Nm."""

        from lcapy.oneport import V, Z

        Voc = self.Voc(Np, Nm).laplace()
        Zoc = self.impedance(Np, Nm)

        return V(Voc) + Z(Zoc)

    def norton(self, Np, Nm):
        """Return s-domain Norton model between nodes Np and Nm."""

        from lcapy.oneport import I, Y

        Isc = self.Isc(Np, Nm).laplace()
        Ysc = self.admittance(Np, Nm)
        
        return I(Isc) | Y(Ysc)

    def admittance(self, Np, Nm):
        """Return s-domain admittance between nodes Np and Nm with independent 
        sources killed.

        """

        new = self.kill()

        # Connect 1 V s-domain voltage source between nodes and
        # measure current.
        new._add('Vin_ %d %d {DiracDelta(t)}' % (Np, Nm))
        If = -new.Vin_.I
        new.remove('Vin_')

        return Ys(If.s)

    def impedance(self, Np, Nm):
        """Return s-domain impedance between nodes Np and Nm with independent
        sources killed.

        """

        new = self.kill()

        # Connect 1 A s-domain current source between nodes and
        # measure voltage.
        new._add('Iin_ %d %d {DiracDelta(t)}' % (Np, Nm))
        Vf = new.Voc(Np, Nm)
        new.remove('Iin_')

        return Zs(Vf.s)

    def transfer(self, N1p, N1m, N2p, N2m):
        """Create voltage transfer function V2 / V1 where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed."""

        new = self.kill()
        new._add('V1_ %d %d {DiracDelta(t)}' % (N1p, N1m))

        V2 = new.Voc(N2p, N2m)
        V1 = new.V1_.V

        return Hs(V2.s / V1.s, causal=True)

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
            self.add('V1_ %d %d {DiracDelta(t)}' % (N1p, N1m))

            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = Hs(self.V1_.V.s / self.Voc(N2p, N2m).s)

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = Zs(self.V1_.V.s / self.Isc(N2p, N2m).s)

            self.remove('V1_')

            self.add('I1_ %d %d {DiracDelta(t)}' % (N1p, N1m))

            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure I2 with port 2 open-circuit
            A21 = Ys(-self.I['I1_'].s / self.Voc(N2p, N2m).s)

            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure I2 with port 2 short-circuit
            A22 = Hs(-self.I['I1_'].s / self.Isc(N2p, N2m).s)

            self.remove('I1_')
            return AMatrix(A11, A12, A21, A22)

        except ValueError:
            raise ValueError('Cannot create A matrix')

    def select(self, sourcenames, kind):
        """Return new netlist with transform domain kind selected for
        specified source.  Sources not in sourcenames are set to zero."""

        new = self._new()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            if cpt.name in sourcenames:
                net = cpt.select(kind)
            elif cpt.source:
                net = cpt.zero()
            elif kind != 'ivp':
                net = cpt.kill_initial()
            else:
                net = str(cpt)
            new._add(net)
        return new        

    def _kill(self, sourcenames):

        new = self._new()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            if cpt.name in self.control_sources:
                net = cpt.zero()                
            elif cpt.name in sourcenames:
                net = cpt.kill()
            elif 'ICs' in sourcenames:
                net = cpt.kill_initial()
            else:
                net = str(cpt)
            new._add(net)
        return new        

    def kill_except(self, *args):
        """Return a new circuit with all but the specified sources killed;
        i.e., make the voltage sources short-circuits and the current
        sources open-circuits.  If no sources are specified, all
        independent sources (including initial conditions) are killed.

        """

        for arg in args:
            if arg not in self.independent_sources and arg != 'ICs':
                raise ValueError('Element %s is not a known source' % arg)
        sources = []
        for source in self.independent_sources:
            if source not in args:
                sources.append(source)
        if 'ICs' not in args:
            sources.append('ICs')            
            
        return self._kill(sources)

    def kill(self, *args):
        """Return a new circuit with the specified sources killed; i.e., make
        the voltage sources short-circuits and the current sources
        open-circuits.  To kill initial conditions, specify `ICs'.  If
        no sources are specified, all independent sources (including
        initial conditions) are killed.

        """

        if len(args) == 0:
            return self.kill_except()

        sources = []
        for arg in args:
            if arg == 'ICs':
                sources.append(arg)
            elif arg in self.independent_sources:
                sources.append(arg)
            else:
                raise ValueError('Element %s is not a known source' % arg)

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
        """"Create s-domain model."""

        new = self._new()
        new.opts = copy(self.opts)

        for cpt in self._elements.values():
            net = cpt.s_model(var)
            new._add(net)
        return new        

    def ac_model(self, var=omega):
        """"Create AC model for specified angular frequency (default
        omega)."""
        
        return self.s_model(j * var)

    def noise_model(self):
        """"Create noise model where resistors are converted into a series
        combination of an ideal resistor and a noise voltage
        source."""
        
        return self.noisy()
    
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
    def has_s(self):
        """Return True if any independent source has an s-domain component."""

        independent_sources = self.independent_sources
        for cpt in independent_sources.values():
            if cpt.has_s:
                return True
        return False

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

        # TODO, ignore zero sources
        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.source)

    def independent_source_groups(self):
        """Return dictionary of source groups.  Each group is a collection of
        sources that can be analysed at the same time.  Noise sources
        have separate groups since they are assumed uncorrelated.

        """

        noise_source_count = 0
        
        groups = {}
        for key, elt in self.elements.items():
            if not elt.source:
                continue
            cpt = elt.cpt
            if cpt.voltage_source:
                cpt_kinds = cpt.Voc.keys()
            else:
                cpt_kinds = cpt.Isc.keys()                
            for cpt_kind in cpt_kinds:
                if cpt_kind == 'n':
                    noise_source_count += 1
                    cpt_kind = 'n%d' % noise_source_count
                if cpt_kind not in groups:
                    groups[cpt_kind] = []
                groups[cpt_kind].append(key)

        return groups    

    @property
    def control_sources(self):
        """Return dictionary of voltage sources required to specify control
        current for CCVS and CCCS components."""

        result = {}

        for key, cpt in self.elements.items():
            if cpt.need_control_current:
                result[cpt.args[0]] = self.elements[cpt.args[0]]
        return result


class Transformdomains(dict):

    def __getattr__(self, attr):
        if attr not in self:
            raise AttributeError('Unknown attribute %s' % attr)
        return self[attr]

    def __getitem__(self, key):
        if key == 'w':
            key = omega
        # This allows a[omega] to work if omega used as key
        # instead of 'omega'.
        if hasattr(key, 'expr'):
            key = key.expr
        return super(Transformdomains, self).__getitem__(key)

    
class Netlist(NetlistMixin, NetfileMixin):
    """This classes handles a generic netlist with multiple sources.
    During analysis, subnetlists are created for each source kind (dc,
    ac, etc).  Since linearity is assumed, superposition is
    employed.

    """

    def __init__(self, filename=None, context=None):

        super (Netlist, self).__init__(filename, context)
        self._invalidate()
        self.kind = 'super'

    def _invalidate(self):

        for attr in ('_sch', '_sub', '_Vdict', '_Idict'):
            try:
                delattr(self, attr)
            except:
                pass

    @property
    def sub(self):
        """Return dictionary of subnetlists keyed by transform domain kind.
        Note, the subnetlists are not created until a specific one is
        selected.

        """

        if hasattr(self, '_sub'):
            return self._sub

        groups = self.independent_source_groups()        

        if self.is_ivp:

            def namelist(elements):
                return ', '.join([elt for elt in elements])

            if self.missing_ic != {}:
                print('Warning: missing initial conditions for %s' %
                      namelist(self.missing_ic))

            newgroups = {'ivp' : []}
            for key, sources in groups.items():
                if isinstance(key, str) and key[0] == 'n':
                    print('Warning: ignoring noise source %'
                          ' for initial value problem' % sources)
                else:
                    newgroups['ivp'] += sources
            groups = newgroups
        
        self._sub = Transformdomains()

        for key, sources in groups.items():
            kind = key
            if isinstance(kind, str) and kind[0] == 'n':
                kind = 'n'
            self._sub[key] = SubNetlist(self, sources, kind)

        return self._sub

    @property
    def kinds(self):
        """Return list of transform domain kinds."""
        return list(self.sub.keys())
    
    @property
    def Vdict(self):
        """Return dictionary of node voltages for each transform domain"""

        try:
            return self._Vdict
        except AttributeError:
            pass        

        result = Nodedict()
        for sub in self.sub.values():
            for node, value in sub.Vdict.items():
                if node not in result:
                    result[node] = Vsuper()
                result[node].add(value)

        self._Vdict = result
        return result

    @property
    def Idict(self):
        """Return dictionary of branch currents for each transform domain"""

        try:
            return self._Idict
        except AttributeError:        
            pass

        result = Branchdict()
        for sub in self.sub.values():
            for node, value in sub.Idict.items():
                if node not in result:
                    result[node] = Isuper()
                result[node].add(value)
        self._Idict = result                    
        return result    

    def get_I(self, name):
        """Current through component"""

        result = Vsuper()
        for sub in self.sub.values():
            I = sub.get_I(name)
            result.add(I)
        result = result.canonical()            
        return result

    def get_i(self, name):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def get_Vd(self, Np, Nm):
        """Voltage drop between nodes"""

        if isinstance(Nm, int):
            Nm = '%s' % Nm
        if isinstance(Np, int):
            Np = '%s' % Np            

        result = Vsuper()
        for sub in self.sub.values():
            Vd = sub.get_Vd(Np, Nm)
            result.add(Vd)
        result = result.canonical()
        return result

    def get_vd(self, Np, Nm):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(Np, Nm).time()

    
class SubNetlist(NetlistMixin, MNA):

    def __new__(cls, netlist, sources, kind):

        # The unwanted sources are zeroed so that we can still refer
        # to them by name, say when wanting the current through them.
        obj = netlist.select(sources, kind=kind)
        obj.kind = kind
        obj.__class__ = cls
        return obj

    def __init__(cls, netlist, sources, kind):
        pass

    def get_I(self, name):
        """Current through component"""

        self._solve()
        return self._Idict[name].canonical()

    def get_i(self):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def get_Vd(self, Np, Nm):
        """Voltage drop between nodes"""

        self._solve()
        return (self._Vdict[Np] - self._Vdict[Nm]).canonical()

    def get_vd(self, Np, Nm):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(Np, Nm).time()

"""
This module provides circuit analysis using modified nodal analysis
(MNA).

The circuit is described using netlists, similar to SPICE, but with
arbitrary node names (except for the ground node which is labelled 0).
The netlists can be loaded from a file or created at run-time.  For
example:

>>> from lcapy import Circuit
>>> cct = Circuit('Voltage divider')
>>> cct.add('V_s fred 0')
>>> cct.add('R_a fred 1')
>>> cct.add('R_b 1 0')
>>> cct.V.pprint()
>>> cct.I.pprint()

cct.V is a directory of the nodal voltages keyed by the node names.
If the nodes are not integers, they need to specified as strings.
cct.I is a directory of the currents through the components keyed by
the component names.  For example,

>>> cct.V['fred'].pprint()
>>> cct.V[1].pprint()
>>> cct.I['R1'].pprint()


Copyright 2014, 2015 Michael Hayes, UCECE
"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from lcapy.core import pprint, Hs, Zs, Ys, Expr, tsym, s, j, omega
from lcapy.oneport import V, I, v, i, Vdc, Idc, Vac, Iac, Vstep, Istep
from lcapy.oneport import Vacstep, Iacstep
from lcapy.oneport import R, L, C, G, Y, Z
from lcapy.twoport import AMatrix, TwoPortBModel
from schematic import Schematic, Opts
from mna import MNA, VCVS, TF, K, TP, Dummy
import re
from copy import copy


__all__ = ('Circuit', )


cpt_types = ['C',  # Capacitor
             'D',  # Diode (not supported)
             'E',  # VCVS
             'F',  # CCCS (not supported yet, can be handled by G)
             'G',  # VCCS (not supported yet)
             'H',  # CCVS (not supported yet, can be handled by E)
             'I',  # Current
             'J',  # JFET (not supported)
             'K',  # Mutual inductance
             'L',  # Inductor
             'M',  # MOSFET (not supported)
             'O',  # Open-circuit
             'P',  # Port (open-circuit)
             'Q',  # BJT (not supported)
             'R',  # Resistor
             'TF',  # Ideal transformer (even works at DC!)
             'TP',  # Two-port (not supported yet)
             'V',  # Voltage
             'W',  # Wire (short-circuit)
             'Y',  # Admittance
             'Z',  # Impedance
             ]

# Note, SPICE netlists are usually case insensitive
# Perhaps prefix mechanical components with M?  But this will make
# crappy component labels.
mech_cpt_types = ['d',  # Dashpot (damper, resistance)  perhaps b?
                  'f',  # Force source
                  'k',  # Spring
                  'm',  # Mass
                  'u',  # Velocity source
                  ]

class Ldict(dict):

    """Lazy dictionary for inverse Laplace"""

    def __init__(self, Vdict):

        super(Ldict, self).__init__()

        self.Vdict = Vdict


    def __getitem__(self, key):

        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key
        
        # Note, need to use keys method to catch branch names.
        if (key not in self) and (key in self.Vdict.keys()):
            v = self.Vdict[key].inverse_laplace()
            self[key] = v
            return v

        return super(Ldict, self).__getitem__(key)

    def keys(self):

        return self.Vdict.keys()


cpt_type_map = {'R': R, 'C': C, 'L': L, 'Z': Z, 'Y': Y,
                'Vac': Vac, 'Vdc': Vdc,
                'Iac': Iac, 'Idc': Idc,
                'Vacstep': Vacstep, 'Vstep': Vstep,
                'Iacstep': Iacstep, 'Istep': Istep,
                'Vimpulse': V, 'Iimpulse': I,
                'Vs': V, 'Is': I,
                'V': V, 'I': I, 'v': v, 'i': i,
                'O' : None, 'P': None, 'W': None,
                'E': VCVS, 'TF': TF, 'TP': TP, 'K': K,
                'D' : Dummy, 'J' : Dummy, 'M': Dummy, 'Q': Dummy,
                'opamp': VCVS}


# Regular expression alternate matches stop with first match so need
# to have longer names first.
cpt_types.sort(lambda x, y: cmp(len(y), len(x)))
cpt_type_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(cpt_types))


class Node(object):

    def __init__(self, cct, name):

        self.cct = cct
        self.name = name
        self.pos = None
        self.port = False
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.list = []

    @property
    def V(self):
        """Node voltage with respect to ground"""

        return self.cct.V[self.name]

    @property
    def v(self):
        """Node time-domain voltage with respect to ground"""

        return self.cct.v[self.name]

    def append(self, elt):

        if elt.cpt_type in ('P', ):
            self.port = True

        self.list.append(elt)


class NetElement(object):

    def __init__(self, cct, name, node1, node2, *args, **opts):

        match = cpt_type_pattern.match(name)

        if not match:
            raise ValueError('Unknown schematic component %s' % name)

        # Circuitikz does not like a . in a name
        if node1.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % node1)
        if node2.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % node2)

        cpt_type = match.groups()[0]
        cpt_id = match.groups()[1]

        # Default value is the component name
        value = cpt_type
        if len(cpt_id) > 0:
            value += '_' + cpt_id
        else:
            if cpt_type == 'W':
                # Automatically enumerate wires to avoid conflict.
                cct.wire_counter += 1
                name = cpt_type + '#%d' % cct.wire_counter

        self.cct = cct
        self.name = name
        self.opts = Opts(opts)
        self.nodes = (node1, node2)
        self.args = args

        # Handle special cases for voltage and current sources.
        # Perhaps could generalise for impulse, step, ramp sources
        # although these can be specified symbolically, for example,
        # v1 1 0 t*Heaviside(t)
        # The only gnarly bit is that the expression cannot contain spaces.

        if cpt_type == 'TP' and len(args) != 5:
            raise ValueError('TP component requires 5 args')

        if args != ():
            if cpt_type in ('V', 'I') and args[0] in (
                    'ac', 'dc', 'step', 'acstep', 'impulse', 's'):
                cpt_type = cpt_type + args[0]
                args = args[1:]
            elif cpt_type == 'E' and args[0] == 'opamp':
                cpt_type = 'opamp'
                args = args[1:]

        if cpt_type in ('E', 'F', 'G', 'H', 'TF', 'TP', 'opamp'):
            if len(args) < 2:
                raise ValueError(
                    'Component type %s requires 4 nodes' % cpt_type)
            self.nodes += (args[0], args[1])
            args = args[2:]

        self.cpt_type = cpt_type

        if args != () and cpt_type in ('V', 'I'):
            # If have a t-domain expression, use v and i.
            expr = Expr(args[0], cache=False)
            if expr.expr.find(tsym) != set():
                cpt_type = 'v' if cpt_type == 'V' else 'i'

        try:
            foo = cpt_type_map[cpt_type]
            # Ignore ports and wires, etc.
            if foo is None:
                self.cpt = None
                return

        except KeyError:
            raise ValueError('Unknown component kind for %s' % name)

        if len(args) == 0:
            # Ensure symbol uppercase for s-domain value.
            if cpt_type in ('Vdc', 'Vac', 'Idc', 'Iac'):
                value = value.capitalize()

            args = (value, )

        cpt = foo(*args)
        self.cpt = cpt

    def __repr__(self):

        args = [self.name] + list(self.nodes) + list(self.args)
        str = ', '.join([arg.__str__() for arg in args])
        return 'NetElement(%s)' % str

    def __str__(self):

        args = (self.name, ) + self.nodes[0:2] + self.args
        return ' '.join(['%s' % arg for arg in args])

    @property
    def _is_dummy(self):

        return self.cpt_type in ('O', 'P', 'W')

    @property
    def _is_V(self):

        return isinstance(self.cpt, (V, Vdc, Vac, Vstep, Vacstep, VCVS, TF))

    @property
    def _is_I(self):

        return isinstance(self.cpt, (I, Idc, Iac, Istep, Iacstep))

    @property
    def _is_RC(self):

        return isinstance(self.cpt, (R, G, C))

    @property
    def _is_L(self):

        return isinstance(self.cpt, L)

    @property
    def _is_K(self):

        return isinstance(self.cpt, K)

    @property
    def I(self):
        """Current through element"""

        return self.cct.I[self.name]

    @property
    def i(self):
        """Time-domain current through element"""

        return self.cct.i[self.name]

    @property
    def V(self):
        """Voltage drop across element"""

        return self.cct.V[self.name]

    @property
    def v(self):
        """Time-domain voltage drop across element"""

        return self.cct.v[self.name]

    @property
    def Y(self):
        """Admittance"""
        
        return self.cpt.Y

    @property
    def Z(self):
        """Impedance"""
        
        return self.cpt.Z


class Netlist(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.sources = {}
        self.nodes = {}
        # Shared nodes (with same voltage)
        self.snodes = {}

        self.wire_counter = 0

        # TODO, decouple from Schematic
        self.kwargs = {'draw_nodes': 'primary',
                       'label_values': True,
                       'label_ids': False,
                       'label_nodes': 'primary'}

        self._MNA = None

        if filename is not None:
            self.netfile_add(filename)

    def __getitem__(self, name):
        """Return element or node by name"""

        # If name is an integer, convert to a string.
        if isinstance(name, int):
            name = '%d' % name

        if name in self.nodes:
            return self.nodes[name]

        if name in self.elements:
            return self.elements[name]

        raise ValueError('Unknown element or node name %s' % name)


    def __getattr__(self, attr):
        """Return element or node by name"""

        # This gets called if there is no explicit attribute attr for
        # this instance.  This is primarily for accessing elements
        # and non-numerical node names.

        return self.__getitem__(attr)


    def netfile_add(self, filename):
        """Add the nets from file with specified filename"""

        file = open(filename, 'r')

        lines = file.readlines()

        for line in lines:
            self.add(line)

    def netlist(self, full=False):
        """Return the current netlist"""

        lines = ''
        for key, elt in self.elements.iteritems():
            newelt = copy(elt)

            if not full:
                newelt.nodes = tuple([self.node_map[node]
                                      for node in elt.nodes])
                if elt._is_dummy:
                    continue

            line = newelt.__str__()
            if full:
                optstr = newelt.opts.format()
                if optstr != '':
                    line += ' ; ' + optstr

            lines += line + '\n'

        return lines

    def _node_add(self, node, elt):

        if node not in self.nodes:
            self.nodes[node] = Node(self, node)
        self.nodes[node].append(elt)

        vnode = self.nodes[node].rootname

        if vnode not in self.snodes:
            self.snodes[vnode] = []

        if node not in self.snodes[vnode]:
            self.snodes[vnode].append(node)

    def _invalidate(self):

        self._MNA = None

        for attr in ('_sch', ):
            if hasattr(self, attr):
                delattr(self, attr)

    def _elt_add(self, elt):

        if elt.name in self.elements:
            print('Overriding component %s' % elt.name)
            # Need to search lists and update component.

        self.elements[elt.name] = elt

        if elt._is_I or elt._is_V:
            self.sources[elt.name] = elt

        # Ignore nodes for mutual inductance.
        if elt.cpt_type == 'K':
            return

        self._node_add(elt.nodes[0], elt)
        self._node_add(elt.nodes[1], elt)

    def net_parse(self, string, *args):

        fields = string.split(';')
        string = fields[1].strip() if len(fields) > 1 else ''
        if string != '':
            self.hints = True

        opts = Opts(string)

        parts = tuple(re.split(r'[\s]+', fields[0].strip()))
        elt = NetElement(self, *(parts + args), **opts)
        return elt

    def add(self, string, *args):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        # Ignore comments
        string = string.strip()
        if string == '' or string[0] in ('#', '%'):
            return

        if '\n' in string:
            lines = string.split('\n')
            for line in lines:
                self.add(line)
            return

        self._invalidate()

        if string[0] == ';':
            keypairs = string[1:].split(',')
            for keypair in keypairs:
                fields = keypair.split('=')
                key = fields[0].strip()
                arg = fields[1].strip() if len(fields) > 1 else ''
                if arg.lower() == 'false':
                    arg = False
                elif arg.lower() == 'true':
                    arg = True
                self.kwargs[key] = arg
            return

        elt = self.net_parse(string, *args)

        self._elt_add(elt)

    def net_add(self, string, *args):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.

        Note, this method has been superseded by add.
        """

        self.add(string, *args)

    def remove(self, name):
        """Remove specified element"""

        self._invalidate()

        if name not in self.elements:
            raise ValueError('Unknown component: ' + name)
        self.elements.pop(name)

    def _make_node(self):
        """Create a dummy node"""

        if not hasattr(self, '_node_counter'):
            self._node_counter = 0
        self._node_counter += 1
        return '_%d' % self._node_counter

    def _make_open(self, node1, node2, opts):
        """Create a dummy open-circuit"""

        if not hasattr(self, '_open_counter'):
            self._open_counter = 0
        self._open_counter += 1

        opts.strip_current_labels()
        opts.strip_labels()

        net = 'P#%d %s %s ; %s' % (
            self._open_counter, node1, node2, opts.format())

        return self.net_parse(net)

    def _make_short(self, node1, node2, opts):
        """Create a dummy short-circuit"""

        if not hasattr(self, '_short_counter'):
            self._short_counter = 0
        self._short_counter += 1

        opts.strip_voltage_labels()
        opts.strip_labels()

        net = 'W#%d %s %s ; %s' % (
            self._short_counter, node1, node2, opts.format())

        return self.net_parse(net)

    def _make_Z(self, name, node1, node2, value, opts):
        """Create a dummy impedance"""

        net = 'Z%s %s %s %s; %s' % (name,
                                    node1, node2, value, opts.format())

        return self.net_parse(net)

    def _make_V(self, node1, node2, value, opts):
        """Create a dummy s-domain voltage source"""

        if not hasattr(self, '_V_counter'):
            self._V_counter = 0
        self._V_counter += 1

        net = 'V#%d %s %s s %s; %s' % (
            self._V_counter, node1, node2, value, opts.format())

        return self.net_parse(net)

    def _make_I(self, node1, node2, value, opts):
        """Create a dummy s-domain current source"""

        if not hasattr(self, '_I_counter'):
            self._I_counter = 0
        self._I_counter += 1

        net = 'I#%d %s %s s %s; %s' % (
            self._I_counter, node1, node2, value, opts.format())

        return self.net_parse(net)

    @property
    def MNA(self):
        """Return results from modified nodal analysis (MNA) of the circuit.

        Note, the voltages and currents are lazily determined when
        requested. """

        if self._MNA is None:

            # TODO: think this out.  When a circuit is converted
            # to a s-domain model we get Z (and perhaps Y) components.
            # We also loose the ability to determine the voltage
            # across a capacitor or inductor since they get split
            # into a Thevenin model and renamed.
            if hasattr(self, '_s_model'):
                raise RuntimeError('Cannot analyse s-domain model')

            self._MNA = MNA(self.elements, self.nodes, self.snodes)
        return self._MNA

    @property
    def V(self):
        """Return dictionary of s-domain node voltages indexed by node name
        and voltage differences indexed by branch name"""

        return self.MNA.V

    @property
    def I(self):
        """Return dictionary of s-domain branch currents
        indexed by component name"""

        return self.MNA.I

    @property
    def v(self):
        """Return dictionary of t-domain node voltages indexed by node name
        and voltage differences indexed by branch name"""

        if not hasattr(self, '_v'):
            self._v = Ldict(self.V)

        return self._v

    @property
    def i(self):
        """Return dictionary of t-domain branch currents indexed
        by component name"""

        if not hasattr(self, '_i'):
            self._i = Ldict(self.I)

        return self._i

    def Voc(self, Np, Nm):
        """Return open-circuit s-domain voltage between nodes Np and Nm."""

        return self.V[Np] - self.V[Nm]

    def voc(self, Np, Nm):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).inverse_laplace()

    def Isc(self, Np, Nm):
        """Return short-circuit s-domain current between nodes Np and Nm."""

        self.add('Vshort_ %d %d' % (Np, Nm), 0)

        Isc = self.I['Vshort_']
        self.remove('Vshort_')

        return Isc

    def isc(self, Np, Nm):
        """Return short-circuit t-domain current between nodes Np and Nm."""

        return self.Isc(Np, Nm).inverse_laplace()

    def thevenin(self, Np, Nm):
        """Return Thevenin model between nodes Np and Nm"""

        Voc = self.Voc(Np, Nm)

        return V(Voc) + Z(self.impedance(Np, Nm))

    def norton(self, Np, Nm):
        """Return Norton model between nodes Np and Nm"""

        Isc = self.Isc(Np, Nm)

        return I(Isc) | Y(self.admittance(Np, Nm))

    def admittance(self, Np, Nm):
        """Return admittance between nodes Np and Nm with sources killed.

        """

        new = self.kill()

        # Connect 1 V s-domain voltage source between nodes and
        # measure current.
        new.add('Vin_ %d %d s 1' % (Np, Nm))
        If = -new.I['Vin_']
        new.remove('Vin_')

        return Ys(If)

    def impedance(self, Np, Nm):
        """Return impedance between nodes Np and Nm with sources killed.

        """

        new = self.kill()

        # Connect 1 A s-domain current source between nodes and
        # measure voltage.
        new.add('Iin_ %d %d s 1' % (Np, Nm))
        Vf = new.Voc(Np, Nm)
        new.remove('Iin_')

        return Zs(Vf)

    def Y(self, Np, Nm):
        """Return admittance between nodes Np and Nm with sources killed.

        """

        return self.admittance(Np, Nm)

    def Z(self, Np, Nm):
        """Return impedance between nodes Np and Nm with sources killed.

        """

        return self.impedance(Np, Nm)

    def transfer(self, N1p, N1m, N2p, N2m):
        """Create voltage transfer function V2 / V1 where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, sources are killed."""

        new = self.kill()
        new.add('V1_ %d %d impulse' % (N1p, N1m))

        H = Hs(new.Voc(N2p, N2m) / new.V['V1_'])

        return H

    def Amatrix(self, N1p, N1m, N2p, N2m):
        """Create A matrix from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """

        if self.Voc(N1p, N1m) != 0 or self.Voc(N2p, N2m) != 0:
            raise ValueError('Network contains independent sources')

        try:

            self.add('V1_ %d %d impulse' % (N1p, N1m))

            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = Hs(self.V['V1_'] / self.Voc(N2p, N2m))

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = Zs(self.V['V1_'] / self.Isc(N2p, N2m))

            self.remove('V1_')

            self.add('I1_ %d %d impulse' % (N1p, N1m))

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

        new = Circuit()

        for key, elt in self.elements.iteritems():
            if key in sourcenames:
                if elt._is_I:
                    newelt = self._make_open(elt.nodes[0], elt.nodes[1], elt.opts)
                elif elt._is_V:
                    newelt = self._make_short(elt.nodes[0], elt.nodes[1], elt.opts)
                else:
                    raise ValueError('Element %s is not a source' % key)
            else:
                newelt = copy(elt)             
                newelt.cct = new
   
            new._elt_add(newelt)

        return new        

    def kill_except(self, *args):
        """Return a new circuit with all but the specified sources killed;
        i.e., make the voltage sources short-circuits and the current
        sources open-circuits.  If no sources are specified, all are
        killed.

        """

        for arg in args:
            if arg not in self.sources:
                raise ValueError('Element %s is not a known source' % arg)
        sources = []
        for key, source in self.sources.iteritems():
            if key not in args:
                sources.append(key)
        return self._kill(sources)

    def kill(self, *args):
        """Return a new circuit with the specified sources killed; i.e., make
        the voltage sources short-circuits and the current sources
        open-circuits.  If no sources are specified, all are killed.

        """

        if len(args) == 0:
            return self.kill_except()

        sources = []
        for arg in args:
            if arg not in self.sources:
                raise ValueError('Element %s is not a known source' % arg)
            sources.append(self.sources[arg].name)

        return self._kill(sources)

    def twoport(self, N1p, N1m, N2p, N2m):
        """Create twoport model from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """

        V2b = self.Voc(N2p, N2m)
        I2b = self.Isc(N2p, N2m)

        A = self.kill().Amatrix(N1p, N1m, N2p, N2m)

        return TwoPortBModel(A.B, V2b, I2b)

    @property
    def sch(self):

        if hasattr(self, '_sch'):
            return self._sch

        netlist = self.netlist(full=True)

        sch = Schematic()

        for net in netlist.split('\n')[0:-1]:
            sch.add(net)

        self._sch = sch
        return sch

    def pre_initial_model(self):
        """Generate circuit model for determining the pre-initial
        conditions."""

        new_cct = self.__class__()

        for key, elt in self.elements.iteritems():

            # Assume initial C voltage and L current is zero.

            if elt.cpt_type in ('V', 'I', 'Vac', 'Iac'):
                print('Cannot determine pre-initial condition for %s'
                      ', assuming 0' % elt.name)

            # v and i should be evaluated to determine the value at 0 - eps.
            if elt.cpt_type in ('v', 'i'):
                print('Cannot determine pre-initial condition for %s'
                      ', assuming 0' % elt.name)

            if elt.cpt_type in ('C', 'Istep', 'Iacstep', 'I', 'i',
                                'Iac', 'Iimpulse'):
                elt = self._make_open(elt.nodes[0], elt.nodes[1], elt.opts)
            elif elt.cpt_type in ('L', 'Vstep', 'Vacstep', 'V', 'v',
                                  'Vac', 'Vimpulse'):
                elt = self._make_short(elt.nodes[0], elt.nodes[1], elt.opts)
            new_cct._elt_add(elt)

        return new_cct

    def _model(self, var=None):

        cct = Circuit()
        cct._s_model = True

        for key, elt in self.elements.iteritems():

            new_elt = copy(elt)

            cpt_type = elt.cpt_type

            if cpt_type in ('C', 'L', 'R'):
                new_elt = self._make_Z(elt.name, elt.nodes[0], elt.nodes[1],
                                       elt.cpt.Z(var), elt.opts)
            elif cpt_type in ('V', 'Vdc', 'Vac', 'Vimpulse',
                              'Vstep', 'Vacstep'):
                new_elt = self._make_V(
                    elt.nodes[0], elt.nodes[1], elt.cpt.V(var), elt.opts)
            elif cpt_type in ('I', 'Idc', 'Iac', 'Iimpulse',
                              'Istep', 'Iacstep'):
                new_elt = self._make_I(
                    elt.nodes[0], elt.nodes[1], elt.cpt.I(var), elt.opts)

            if cpt_type in ('C', 'L', 'R') and elt.cpt.V != 0:

                dummy_node = self._make_node()

                velt = self._make_V(
                    dummy_node, elt.nodes[1], elt.cpt.V(var), elt.opts)
                new_elt.nodes = (elt.nodes[0], dummy_node)

                # Strip voltage labels. 
                voltage_opts = new_elt.opts.strip_voltage_labels()

                # Strip voltage and current labels from voltage source.
                velt.opts.strip_all_labels()

                cct._elt_add(velt)

                if voltage_opts != {}:
                    opts = elt.opts.copy()
                    opts.strip_current_labels()
                    # Need to convert voltage labels to s-domain.
                    # v(t) -> V(s)
                    # v_C -> V_C
                    # v_L(t) -> V_L(s)
                    for opt, val in voltage_opts.iteritems():
                        opts[opt] = val.capitalize()

                    open_elt = self._make_open(elt.nodes[0], elt.nodes[1], opts)
                    cct._elt_add(open_elt)

            # Make voltage and current labels uppercase.
            for opt, val in new_elt.opts.strip_voltage_labels().iteritems():
                new_elt.opts[opt] = val.capitalize()            
            for opt, val in new_elt.opts.strip_current_labels().iteritems():
                new_elt.opts[opt] = val.capitalize()            
            cct._elt_add(new_elt)

        return cct

    def ac_model(self):
        return self._model(j * omega)


    def s_model(self):
        return self._model(s)


    def draw(self, filename=None, label_values=None, draw_nodes=None,
             label_nodes=None, label_ids=None,
             s_model=False, args=None, scale=1, stretch=1,
             **kwargs):

        cct = self
        if s_model:
            cct = cct.s_model()

        kwargs2 = copy(self.kwargs)
        if draw_nodes is not None:
            kwargs2['draw_nodes'] = draw_nodes
        if label_nodes is not None:
            kwargs2['label_nodes'] = label_nodes
        if label_values is not None:
            kwargs2['label_values'] = label_values
        if label_ids is not None:
            kwargs2['label_ids'] = label_ids

        for key, arg in kwargs.iteritems():
            kwargs2[key] = arg

        return cct.sch.draw(filename=filename, args=args,
                            scale=scale, stretch=stretch,
                            **kwargs2)


class Circuit(Netlist):

    def __init__(self, filename=None):

        super(Circuit, self).__init__(filename)


def test():

    cct = Circuit('Test')

    cct.add('V_s fred 0')
    cct.add('R_a fred bert')
    cct.add('R_b bert 0')

    pprint(cct.V)

    pprint(cct.I)

    return cct

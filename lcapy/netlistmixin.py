"""This module provides the NetlistMixin class.  This is used for
Netlist and SubNetlist.

Copyright 2020-2021 Michael Hayes, UCECE

"""

from .expr import expr
from .mnacpts import Cpt
from .impedance import impedance
from .admittance import admittance
from .equipotentialnodes import EquipotentialNodes
from .node import Node
from .state import state
from .schematic import Schematic
from .symbols import j, s, omega
from .attrdict import AttrDict
from .netfile import NetfileMixin
from .statespace import StateSpace
from .voltage import Vname
from .current import Iname, current
from .simulator import Simulator
from .netlistnamespace import NetlistNamespace
from .matrix import Matrix

from . import mnacpts
from collections import OrderedDict


class NetlistMixin(object):

    def __init__(self, filename=None, context=None, allow_anon=False):

        self._elements = OrderedDict()
        self.namespaces = {}
        self.nodes = AttrDict()
        if context is None:
            context = state.new_context()
        
        self.context = context
        self.allow_anon = allow_anon
        self._init_parser(mnacpts, allow_anon=allow_anon)

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

        if name in self.namespaces:
            return self.namespaces[name]
        
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

    def pdb(self):
        """Enter the python debugger."""
        
        import pdb; pdb.set_trace()
        return self

    @property
    def symbols(self):
        """Return dictionary of symbols defined in the circuit."""
        
        return self.context.symbols

    @property
    def all_symbols(self):
        """Return dictionary of symbols defined in the circuit and global
        symbols."""

        symbols = self.symbols
        symbols.update(state.global_context.symbols)
        return symbols
    
    @property
    def elements(self):

        if hasattr(self, '_add_elements'):
            if self._elements == {}:
                self._add_elements()

        return self._elements

    def netlist(self):
        """Return the current netlist."""

        return '\n'.join([str(cpt) for cpt in self._elements.values()])

    def has(self, cpt):
        """Return True if cpt in elements."""

        return cpt in self.elements
    
    @property
    def is_connected(self):
        """Return True if all components are connected."""

        return self.cg.is_connected

    def _node_add(self, node, cpt):

        if node not in self.nodes:
            self.nodes[node] = Node(self, node)
        self.nodes[node].append(cpt)

    def _cpt_add(self, cpt):

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

        for node in cpt.nodenames:
            self._node_add(node, cpt)

        self._namespace_add(cpt.namespace)

    def _namespace_add(self, namespace):

        namespace = namespace.strip('.')
        if namespace == '':
            return

        parts = namespace.split('.')
        namespaces = self.namespaces
        namespace = ''
        for part in parts:
            if part not in namespaces:
                namespace += part
                namespaces[part] = NetlistNamespace(namespace, self)
                namespace += '.'
                namespaces = namespaces[part].namespaces
            
    def copy(self):
        """Create a copy of the netlist"""

        new = self._new()

        for cpt in self._elements.values():
            new._add(cpt._copy())
        return new        

    def _new(self):

        from .circuit import Circuit
        from .netlist import Netlist
        
        # TODO.  Copy or share?
        context = self.context
        if self.__class__ == Circuit:
            return Circuit(context=context)
        # If have OnePort, Network, etc., treat as Netlist
        return Netlist(context=context)

    def remove(self, name):
        """Remove specified element or elements specified in list."""

        if isinstance(name, (list, tuple)):
            for name1 in name:
                self.remove(name1)
            return self

        self._invalidate()

        if name not in self._elements:
            raise ValueError('Unknown component: ' + name)
        self._elements.pop(name)
        # TODO, remove nodes that are only connected
        # to this component.
        return self

    @property
    def super_nodes(self):
        """Super nodes are nodes linked by voltage sources."""

        snodes = []

        for elt in self.elements.values():
            if not elt.is_voltage_source:
                continue
            snodes.append(elt.nodenames[0:2])

        from .utils import merge_common

        # Merge supernodes to create super-supernodes...
        return list(merge_common(snodes))

    @property
    def equipotential_nodes(self):
        """Determine nodes connected by wires that are of the same potential.
        This returns a dictionary keyed by the unique node names with
        values being lists of nodes of the same potential."""

        enodes = EquipotentialNodes()
        enodes.add(self.nodes.keys())

        # Then augment with nodes connected by wires.
        for m, elt in enumerate(self.elements.values()):
            if elt.type == 'W':
                enodes.add_wire(*elt.nodenames)
            else:
                for connections in elt.equipotential_nodes:
                    enodes.add_wires([elt.name + '.' + n for n in connections])

        # Alter keys to avoid underscore and to ensure that have a '0'
        # key if possible.
        enodes2 = {}
        for key, nodes in enodes.items():
            nodes = sorted(nodes)
            if '0' in nodes:
                newkey = '0'
            else:
                newkey = nodes[0]
            enodes2[newkey] = nodes
                
        return enodes2

    @property
    def node_map(self):
        """Create dictionary mapping node names to the unique
        equipotential node names."""

        if hasattr(self, '_node_map'):
            return self._node_map

        enodes = self.equipotential_nodes

        # Create inverted dictionary that maps the node names
        # to the equipotential node names.
        node_map = {}        
        for key, nodes in enodes.items():
            for node in nodes:
                node_map[node] = key

        self._node_map = node_map
        return node_map

    def annotate_current(self, cpts, var=None, flow=False, pos=''):
        """Annotate specified list of component names `cpts` with current (or
        flow).  `pos` specifies where to position the labels."""

        new = self._new()        
        for cpt in self._elements.values():
            net = cpt._copy()
            if cpt.name in cpts:
                I = cpt.I
                if var is not None:
                    I = I(var)

                net += ', ' if ';' in net else '; '
                if flow:
                    net += ', f%s=$%s$' % (pos, I.latex())
                else:
                    net += ', i%s=$%s$' % (pos, I.latex())
            new.add(net)
        return new

    def annotate_voltage(self, cpts, var=None, pos=''):
        """Annotate specified list of component names `cpts` with voltage.
        `pos` specifies where to position the labels."""

        new = self._new()                
        for cpt in self._elements.values():
            net = cpt._copy()
            if cpt.name in cpts:
                V = cpt.V
                if var is not None:
                    V = V(var)
                net += ', ' if ';' in net else '; '                    
                net += 'v%s=$%s$' % (pos, V.latex())
            new.add(net)
        return new                
    
    def augment_node_map(self, node_map=None):

        """Create a mapping dict for all nodes."""

        if node_map is None:
            node_map = {}
        
        # It would be desirable to renumber the nodes say from left to
        # right and top to bottom.  The schematic drawing algorithms
        # could help with this since they figure out the node
        # placement. 
        
        enodes = self.equipotential_nodes

        if '0' in self.nodes and '0' not in node_map:
            node_map['0'] = '0'
        
        numbers = []
        for m in range(len(enodes)):
            numbers.append('%s' % (m + 1))
        
        for old, new in node_map.items():
            if old not in self.nodes:
                raise ValueError('Unknown node %s' % old)
            if new in numbers:
                numbers.remove(new)

        # Ensure that enode keys are the nodes to be renamed.
        enodes2 = {}
        for key, nodes in enodes.items():
            newkey = None
            for node in nodes:
                if node not in node_map:
                    continue
                if newkey is not None:
                    # TODO, FIXME
                    raise ValueError('Cannot rename two nodes of same potential')
                newkey = node_map[node]
            if newkey is None:
                newkey = numbers.pop(0)
            enodes2[key] = (newkey, nodes)

        for key, foo in enodes2.items():

            newkey, nodes = foo
            
            snodes = sorted(nodes.copy())

            root = newkey
            snodes.remove(key)
            node_map[key] = root
                
            for m, enode in enumerate(snodes):
                node_map[enode] = root + '_%d' % (m + 1)

        return node_map
    
    def renumber(self, node_map=None):
        """Renumber nodes using specified node_map.  If node_map not specified
        then a mapping is created."""

        if node_map is None:
            node_map = {}
        
        if len(node_map) != len(self.nodes):
            node_map = self.augment_node_map(node_map)

        new = self._new()

        for cpt in self._elements.values():
            new._add(cpt._rename_nodes(node_map))
        return new                

    @property
    def node_list(self):
        """Determine list of sorted unique node names, e.g.,
        ['0', '1', '2']."""

        if hasattr(self, '_node_list'):
            return self._node_list

        # Extract unique nodes.
        node_list = list(self.equipotential_nodes.keys())
        node_list = sorted(node_list)
        # Ensure node '0' is first in the list.
        if '0' in node_list:
            node_list.insert(0, node_list.pop(node_list.index('0')))

        self._node_list = node_list
        return node_list

    @property
    def branch_list(self):
        """Determine list of names of branch elements, e.g.,
        ['C1', 'V1', 'R1', 'R2'b]."""

        if hasattr(self, '_branch_list'):
            return self._branch_list

        self._branch_list = []
        for key, elt in self.elements.items():
            if elt.type not in ('W', 'O', 'P', 'K'):
                self._branch_list.append(elt.name)                
        return self._branch_list

    def _check_nodes(self, *nodes):

        str_nodes = []
        for node in nodes:
            if isinstance(node, int):
                node = '%s' % node
            if node not in self.nodes:
                raise ValueError('Unknown node %s' % node)
            str_nodes.append(node)
        return str_nodes
    
    def _parse_node_args(self, Np, Nm=None):
        
        if Nm is None:
            cpt = self[Np]
            if isinstance(cpt, Node):
                Np, Nm = cpt.name, 0
            else:
                Np, Nm = cpt.nodenames[0:2]
        return Np, Nm
    
    def Voc(self, Np, Nm=None):
        """Return open-circuit transform-domain voltage between nodes Np and
        Nm."""

        return self.get_Vd(Np, Nm)

    def voc(self, Np, Nm=None):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).time()

    def Isc(self, Np, Nm=None):
        """Return short-circuit transform-domain current between nodes Np and
        Nm."""

        Np, Nm = self._parse_node_args(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        
        new = self.copy()
        if new.is_causal:
            new.add('Vshort_ %s %s step 0' % (Np, Nm))
        else:
            new.add('Vshort_ %s %s 0' % (Np, Nm))            

        # Negate current since Vshort is a considered a source.
        Isc = -new.Vshort_.I
        new.remove('Vshort_')
        
        return Isc

    def isc(self, Np, Nm=None):
        """Return short-circuit t-domain current between nodes Np and Nm."""

        return self.Isc(Np, Nm).time()

    def oneport(self, Np, Nm=None):
        """Return oneport object between nodes Np and Nm.  This might be a
        Thevenin network, a Norton network, or a single component.

        If Np is a component name, create model using the component nodes."""

        try:
            return self.norton(Np, Nm).simplify()
        except:
            return self.thevenin(Np, Nm).simplify()
    
    def thevenin(self, Np, Nm=None):
        """Return s-domain Thevenin oneport model between nodes Np and Nm.

        If Np is a component name, create model using the component nodes."""

        from .oneport import V, Z

        Np, Nm = self._parse_node_args(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)        
        Voc = self.Voc(Np, Nm)
        Zoc = self.impedance(Np, Nm)

        # Convert to time-domain to handle arbitrary sources.  Either
        # this or define a way to represent a superposition in a
        # netlist.
        return V(Voc.time()) + Z(Zoc)

    def norton(self, Np, Nm=None):
        """Return s-domain Norton model between nodes Np and Nm.

        If Np is a component name, create model using the component nodes."""

        from .oneport import I, Y

        Np, Nm = self._parse_node_args(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)                        
        Isc = self.Isc(Np, Nm)
        Ysc = self.admittance(Np, Nm)

        # Convert to time-domain to handle arbitrary sources.  Either
        # this or define a way to represent a superposition in a
        # netlist.        
        return I(Isc.time()) | Y(Ysc)

    def match(self, pattern):
        """Return list of components names matching regular 
        expression pattern."""

        from re import match

        elts = []
        for cpt in self._elements.values():
            if match(pattern, cpt.name):
                elts.append(cpt.name)
        return elts

    
    def admittance(self, Np, Nm=None):
        """Return driving-point admittance between nodes
        Np and Nm with independent sources killed and initial
        conditions ignored.  Since the result is causal, the frequency
        domain admittance can be found by substituting j * omega for
        s."""        

        Np, Nm = self._parse_node_args(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        
        new = self.kill()
        if '0' not in new.nodes:
            new.add('W %s 0' % Nm)        

        # Connect 1 V s-domain voltage source between nodes and
        # measure current.
        new._add('Vin_ %s %s {DiracDelta(t)}' % (Np, Nm))
        If = new.Vin_.I
        new.remove('Vin_')

        return admittance(If.laplace().expr)

    def impedance(self, Np, Nm=None):
        """Return driving-point impedance between nodes
        Np and Nm with independent sources killed and initial
        conditions ignored.  Since the result is causal, the frequency
        domain impedance can be found by substituting j * omega for
        s."""

        Np, Nm = self._parse_node_args(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)        
        
        new = self.kill()
        if '0' not in new.nodes:
            new.add('W %s 0' % Nm)

        # Connect 1 A s-domain current source between nodes and
        # measure voltage.
        new._add('Iin_ %s %s {DiracDelta(t)}' % (Np, Nm))
        Vf = new.Voc(Np, Nm)
        new.remove('Iin_')

        return impedance(Vf.laplace().expr)

    def resistance(self, Np, Nm=None):
        """Return resistance between nodes Np and Nm with independent
        sources killed.  The result is in the AC (omega) domain.
        See also conductance, reactance, susceptance."""
        return self.impedance(Np, Nm).R

    def reactance(self, Np, Nm=None):
        """Return reactance between nodes Np and Nm with independent
        sources killed.  The result is in the AC (omega) domain.
        See also conductance, resistance, susceptance."""
        return self.impedance(Np, Nm).X

    def conductance(self, Np, Nm=None):
        """Return conductance (inverse resistance) between nodes Np and Nm
        with independent sources killed.  The result is in the AC (omega)
        domain.    See also resistance, reactance, susceptance."""
        return self.impedance(Np, Nm).G

    def susceptance(self, Np, Nm=None):
        """Return susceptance (inverse reactance) between nodes Np and Nm with
        independent sources killed.  The result is in the AC (omega)
        domain.   See also conductance, reactance, resistance."""
        return self.impedance(Np, Nm).B
        
    def transfer(self, N1p, N1m, N2p=None, N2m=None):
        """Create s-domain voltage transfer function V2(s) / V1(s) where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed and initial conditions
        are ignored.  Since the result is causal, the frequency response
        can be found by substituting j * omega for s.

        An alternative syntax is transfer(cpt1, cpt2)."""

        if N2p is None and N2m is None:
            try:
                N1p = self.elements[N1p]
            except:
                pass
            try:
                N1m = self.elements[N1m]
            except:
                pass            
            
            if N1p not in self.elements.values():
                raise ValueError('Unknown component %s' % N1p)
            if N1m not in self.elements.values():
                raise ValueError('Unknown component %s' % N1m)            
            N2p, N2m = [n.name for n in N1m.nodes[0:2]]
            N1p, N1m = [n.name for n in N1p.nodes[0:2]]
            
        elif N2p is None or N2m is None:
            raise ValueError('Expecting transfer(cpt1, cpt2) or transfer(N1p, N1m, N2p, N2m)')

        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)
        
        new = self.kill()
        if '0' not in new.nodes:
            new.add('W %s 0' % N1m)
        
        new._add('V1_ %s %s {DiracDelta(t)}' % (N1p, N1m))

        V2 = new.Voc(N2p, N2m)
        V1 = new.V1_.V

        H = V2.laplace() / V1.laplace()
        H.causal = True
        return H

    def Aparams(self, N1p, N1m, N2p, N2m):
        """Create A-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Bparams, Gparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """

        from .twoport import AMatrix

        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)        
        net = self.kill()        
        if '0' not in net.nodes:
            net.add('W %s 0' % N1m)                           

        try:
            net.add('V1_ %s %s {DiracDelta(t)}' % (N1p, N1m))

            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = net.V1_.V(s) / net.Voc(N2p, N2m)(s)

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = net.V1_.V(s) / net.Isc(N2p, N2m)(s)

            net.remove('V1_')

            net.add('I1_ %s %s {DiracDelta(t)}' % (N1p, N1m))

            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure V2 with port 2 open-circuit
            try:
                A21 = current(0 * s + 1) / net.Voc(N2p, N2m)(s)
            except ValueError:
                # It is likely there is an open-circuit.                
                net2 = net.copy()
                net2.add('W %s %s' % (N2p, N2m))
                A21 = -net2.I1_.I(s) / net2.Voc(N2p, N2m)(s)
                A21 = 0                

            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure I2 with port 2 short-circuit
            A22 = current(0 * s + 1) / net.Isc(N2p, N2m)(s)

            net.remove('I1_')
            A = AMatrix(((A11, A12), (A21, A22)))
            return A

        except ValueError:
            raise ValueError('Cannot create A matrix')

    def Bparams(self, N1p, N1m, N2p, N2m):
        """Create B-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Aparams, Gparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Bparams

    def Gparams(self, N1p, N1m, N2p, N2m):
        """Create G-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Aparams, Bparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Gparams

    def Hparams(self, N1p, N1m, N2p, N2m):
        """Create H-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Aparams, Bparams, Gparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Hparams

    def Sparams(self, N1p, N1m, N2p, N2m):
        """Create S-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Aparams, Bparams, Gparams, Hparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Sparams

    def Tparams(self, N1p, N1m, N2p, N2m):
        """Create T-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Aparams, Bparams, Gparams, Hparams, Sparams, Yparams, and Zparams.
        """
        return self.Tparams(N1p, N1m, N2p, N2m).Hparams

    def Yparams(self, N1p, N1m, N2p, N2m):
        """Create Y-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Aparams, Bparams, Gparams, Hparams, Sparams, Tparams, and Zparams.
        """
        return self.Zparams(N1p, N1m, N2p, N2m).Yparams            

    def Zparams(self, N1p, N1m, N2p, N2m):
        """Create Z-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also  Aparams, Bparams, Gparams, Hparams, Sparams, Tparams, and Yparams.
        """
        from .twoport import ZMatrix

        # TODO, generalise to multiports.
        
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)        
        net = self.kill()
        if '0' not in net.nodes:
            net.add('W %s 0' % N1m)

        try:
            net.add('I1_ %s %s {DiracDelta(t)}' % (N1p, N1m))

            # Z11 = V1 / I1 with I2 = 0
            # Apply I1 and measure V1 with port 2 open-circuit
            Z11 = impedance(net.Voc(N1p, N1m)(s))

            # Z21 = V2 / I1 with I2 = 0
            # Apply I1 and measure V2 with port 2 open-circuit
            Z21 = impedance(net.Voc(N2p, N2m)(s))          

            net.remove('I1_')

            net.add('I2_ %s %s {DiracDelta(t)}' % (N2p, N2m))

            # Z12 = V1 / I2 with I1 = 0
            # Apply I2 and measure V1 with port 1 open-circuit
            Z12 = impedance(net.Voc(N1p, N1m)(s))

            # Z22 = V2 / I2 with I1 = 0
            # Apply I2 and measure V2 with port 1 open-circuit
            Z22 = impedance(net.Voc(N2p, N2m)(s))          

            net.remove('I2_')            

            Z = ZMatrix(((Z11, Z12), (Z21, Z22)))
            return Z

        except ValueError:
            raise ValueError('Cannot create Z matrix')

    def Yparamsn(self, *nodes):
        """Create Y-parameters for N-port defined by list of node-pairs.

        See also Yparams for a two port.

        """

        nodes = self._check_nodes(*nodes)
        if len(nodes) % 2 == 1:
            raise ValueError('Need an even number of nodes.')
        ports = []
        for m in range(len(nodes) // 2):
            ports.append((nodes[m * 2], nodes[m * 2 + 1]))
        
        net = self.kill()
        if '0' not in net.nodes:
            net.add('W %s 0' % nodes[1])

        try:

            Y = Matrix.zeros(len(ports))
            
            for col in range(len(ports)):

                for row in range(len(ports)):
                    if row == col:
                        net.add('V%d_ %s %s {DiracDelta(t)}' % (row, ports[row][0], ports[row][1]))
                    else:
                        net.add('V%d_ %s %s 0' % (row, ports[row][0], ports[row][1]))                        

                for row in range(len(ports)):
                    Y[row, col] = admittance(net.elements['V%d_' % row].I(s))

                for row in range(len(ports)):                        
                    net.remove('V%d_' % row)
            return Y

        except ValueError:
            raise ValueError('Cannot create Y matrix')

    def Yparams3(self, N1p, N1m, N2p, N2m, N3p, N3m):
        """Create Y-parameters for three-port defined by nodes N1p, N1m, N2p,
        N2m, N3p, and N3m, where:

        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        I3 is the current flowing into N3p and out of N3m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        V3 is V[N3p] - V[N3m]

        See also Yparams for a two port and Yparamsn for an N-port.

        """

        return self.Yparamsn(N1p, N1m, N2p, N2m, N3p, N3m)        

    def Zparamsn(self, *nodes):
        """Create Z-parameters for N-port defined by list of node-pairs.

        See also Zparams for a two port.

        """

        nodes = self._check_nodes(*nodes)
        if len(nodes) % 2 == 1:
            raise ValueError('Need an even number of nodes.')
        ports = []
        for m in range(len(nodes) // 2):
            ports.append((nodes[m * 2], nodes[m * 2 + 1]))
        
        net = self.kill()
        if '0' not in net.nodes:
            net.add('W %s 0' % nodes[1])

        try:

            Z = Matrix.zeros(len(ports))
            
            for col in range(len(ports)):
                net.add('I_ %s %s {DiracDelta(t)}' % (ports[col][0], ports[col][1]))

                for row in range(len(ports)):                
                    Z[row, col] = impedance(net.Voc(ports[row][0], ports[row][1])(s))

                net.remove('I_')
            return Z

        except ValueError:
            raise ValueError('Cannot create Z matrix')        

    def Zparams3(self, N1p, N1m, N2p, N2m, N3p, N3m):
        """Create Z-parameters for three-port defined by nodes N1p, N1m, N2p,
        N2m, N3p, and N3m, where:

        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        I3 is the current flowing into N3p and out of N3m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        V3 is V[N3p] - V[N3m]

        See also Zparams for a two port and Zparamsn for an N-port.

        """

        return self.Zparamsn(N1p, N1m, N2p, N2m, N3p, N3m)

    def save(self, filename):
        """Save netlist to file."""

        f = open(filename, 'w')
        f.writelines(str(self))
        f.close()
        
    def select(self, kind):
        """Return new netlist with transform domain kind selected for
        specified sources in sourcenames. 

        """

        new = self._new()

        for cpt in self._elements.values():
            if cpt.independent_source:
                net = cpt._select(kind)                
            else:
                net = cpt._copy()
            new._add(net)
        return new        

    def _kill(self, sourcenames):

        new = self._new()

        for cpt in self._elements.values():
            if cpt.name in sourcenames:
                if cpt.name in self.control_sources:
                    net = cpt._zero()                
                else:
                    net = cpt._kill()
            else:
                net = cpt._copy()                
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
                raise ValueError('Element %s is not a known independent source' % arg)
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
                raise ValueError('Element %s is not a known independent source' % arg)

        return self._kill(sources)

    def kill_zero(self):
        """Return a new circuit with the independent zero voltage sources and
        zero current sources killed."""

        new = self._new()

        for cpt in self._elements.values():
            if (cpt.independent_source and 
                (cpt.is_voltage_source and cpt.Voc == 0) or
                (cpt.is_current_source and cpt.Isc == 0)):
                net = cpt._kill()
            else:
                net = cpt._copy()
            new._add(net)
        return new            

    def _noisy(self, resistornames, T='T'):

        new = self._new()

        for cpt in self._elements.values():
            if cpt.name in resistornames:
                net = cpt._noisy(T=T)
            else:
                net = cpt._copy()
            new._add(net)
        return new        

    def noisy_except(self, *args, T='T'):
        """Return a new circuit with all but the specified resistors in series
        with noise voltage sources"""

        for arg in args:
            if arg not in self.elements and not self.elements[arg].is_resistor:
                raise ValueError('Element %s is not a known resistor' % arg)
        resistors = []
        for cpt in self.elements.values():
            if cpt.is_resistor and not cpt.is_noiseless and cpt.name not in args:
                resistors.append(cpt.name)
        return self._noisy(resistors, T)

    def noisy(self, *args, T='T'):
        """Return a new circuit with the specified resistors in series
        with noise voltage sources"""

        if len(args) == 0:
            return self.noisy_except(T=T)

        resistors = []
        for arg in args:
            if arg not in self.elements and not self.elements[arg].is_resistor:
                raise ValueError('Element %s is not a known resistor' % arg)
            resistors.append(arg)

        return self._noisy(resistors, T=T)

    @property
    def cg(self):
        """Generate graph for this netlist.   This is cached."""        

        from .circuitgraph import CircuitGraph
        
        if hasattr(self, '_cg'):
            return self._cg

        self._cg = CircuitGraph(self)
        return self._cg
    
    def _potential_combine_names(self):

        names = []
        for name, elt in self.elements.items():
            if elt.type in ('V', 'I', 'R', 'NR', 'C', 'L', 'Y', 'Z'):
                names.append(name)
        return names

    def _find_combine_subsets(self, aset):
        """Return dict of subsets of component names where each subset has the
        same component type."""

        aset = aset.copy()

        subsets = {}
        while aset != set():
            name = aset.pop()
            cpt = self._elements[name]
            aset.add(name)
            
            subset = set()
            for name1 in aset:
                cpt1 = self._elements[name1]                
                if cpt.type == cpt1.type:
                    subset.add(name1)
            aset -= subset
            if len(subset) > 1:
                subsets[cpt.type] = subset
        return subsets

    def _in_parallel_all(self):

        names = self._potential_combine_names()
        lists = []
        while names != []:
            name = names[0]
            parallel = self._in_parallel_names(name)
            if len(parallel) > 1:
                lists.append(parallel)
            for name in parallel:
                try:
                    names.remove(name)
                except:
                    pass
        return lists

    def _in_series_all(self):

        names = self._potential_combine_names()
        lists = []
        while names != []:
            name = names[0]
            series = self._in_series_names(name)
            if len(series) > 1:            
                lists.append(series)
            for name in series:
                try:
                    names.remove(name)
                except:
                    pass                
        return lists    

    def _in_parallel_names(self, cpt=None):

        if isinstance(cpt, Cpt):
            cpt = cpt.name

        return self.cg.in_parallel(cpt)

    def _in_series_names(self, cpt=None):

        if isinstance(cpt, Cpt):
            cpt = cpt.name

        return self.cg.in_series(cpt)
        
    def in_parallel(self, cpt=None):
        """Return set of cpts in parallel with specified cpt.  If no cpt
        specified, return list of sets of parallel cpts."""

        if cpt is None:
            return self._in_parallel_all()
        
        names = self._in_parallel_names(cpt)
        if len(names) < 2:
            return set()
        return names
    
    def in_series(self, cpt=None):
        """Return set of cpts in series with specified cpt.  If no cpt
        specified, return list of sets of series cpts."""

        if cpt is None:
            return self._in_series_all()

        names = self._in_series_names(cpt)
        if len(names) < 2:
            return set()
        return names

    def _do_simplify_combine(self, string, subset, net, explain=False, add=False,
                     series=False):

        if explain:
            print(string % subset)

        subset_list = list(subset)
            
        if add:
            total = expr(0)
            for name in subset_list:
                total += expr(self.elements[name].cpt.args[0])
        else:
            total = expr(0)
            for name in subset_list:
                total += (1 / expr(self.elements[name].cpt.args[0]))
            total = 1 / total

        if explain:
            print('%s combined value = %s' % (subset, total))

        ic = None
        name = subset_list[0]
        elt = self.elements[name]
        if elt.cpt.has_ic:
            ic = expr(0)
            for name1 in subset_list:
                ic += expr(self.elements[name1].cpt.args[1])

            if explain:
                print('%s combined IC = %s' % (subset, ic))                

        newname = self.namer(name[0], 't')                
        net1 = elt._new_value(total, ic)
        parts = net1.split(' ', 1)
        net1 = newname + ' ' + parts[1]
        
        elt = self._parse(net1)
        
        # Overwrite component with one having total value.  _fixup will
        # fix element keys later on once all simplifications are performed.
        net.elements[name] = elt

        for name1 in subset_list[1:]:
            # Replace with wire or open-circuit.
            if series:
                net1 = self.elements[name1]._netmake_W()
            else:
                net1 = self.elements[name1]._netmake_O()
            # Remove component
            elt = self._parse(net1)
            net.elements[name1] = elt
                
        return True

    def _check_ic(self, subset):

        subset = subset.copy()
        name = subset.pop()
        has_ic = self.elements[name].has_ic

        okay = True
        for name1 in subset:
            if self.elements[name1].has_ic != has_ic:
                print('Incompatible initial conditions for %s and %s' % (name, name1))
                okay = False
        if not has_ic:
            return okay
        ic = self.elements[name].cpt.args[1]
        for name1 in subset:
            if self.elements[name1].cpt.args[1] != ic:
                print('Incompatible initial conditions for %s and %s' % (name, name1))
                okay = False

        return okay

    def _fixup(self):
        """Rename keys to fix things up for removed components."""
        
        newelements = OrderedDict()
        for k, v in self.elements.items():
            newelements[v.name] = v
        self._elements = newelements
    
    def _simplify_combine_series(self, cptnames=None, explain=False):

        net = self.copy()
        changed = False

        for aset in net.in_series():
            subsets = net._find_combine_subsets(aset)
            for k, subset in subsets.items():
                if k == 'I':
                    print('Netlist has current sources in series: %s' % subset)
                elif k in ('R', 'NR', 'L', 'V', 'Z'):
                    if k == 'L' and not self._check_ic(subset):
                        continue
                    changed |= self._do_simplify_combine('Can add in series: %s',
                                                         subset, net, explain, True, True)
                elif k in ('C', 'Y'):
                    changed |= self._do_simplify_combine('Can combine in series: %s',
                                                         subset, net, explain, False, True)
                else:
                    raise RuntimeError('Internal error')

        if changed:
            net._fixup()
                
        return net, changed

    def _simplify_combine_parallel(self, cptnames=None, explain=False):

        net = self.copy()
        changed = False        

        for aset in net.in_parallel():
            subsets = net._find_combine_subsets(aset)
            for k, subset in subsets.items():
                if k == 'V':
                    print('Netlist has voltage sources in parallel: %s'% subset)
                elif k in ('R', 'NR', 'L', 'Z'):
                    changed |= self._do_simplify_combine('Can combine in parallel: %s',
                                                         subset, net, explain, False, False)
                elif k in ('C', 'Y', 'I'):
                    if k == 'C' and  not self._check_ic(subset):
                        continue                    
                    changed |= self._do_simplify_combine('Can add in parallel: %s',
                                                         subset, net, explain, True, False)
                else:
                    raise RuntimeError('Internal error')
                
        if changed:
            # TODO, remove dangling wires connected to the removed components.
            net._fixup()
                
        return net, changed

    def _simplify_redundant_series(self, cptnames=None, explain=False):

        net = self.copy()
        changed = False

        for aset in net.in_series():        
            Iname = None
            for name in aset:
                cpt = self._elements[name]                
                if cpt.type == 'I':
                    Iname = name
                    break
            if Iname is not None:
                for name in aset:
                    cpt = self._elements[name]                    
                    if cpt.type != 'I':
                        print('Warning, have redundant %s in series with %s' % (name, Iname))
                
        return net, False        

    def _simplify_redundant_parallel(self, cptnames=None, explain=False):

        net = self.copy()
        changed = False

        for aset in net.in_parallel():
            Vname = None
            for name in aset:
                cpt = self._elements[name]
                if cpt.type == 'V':
                    Vname = name
                    break
            if Vname is not None:
                for name in aset:
                    cpt = self._elements[name]
                    if cpt.type != 'V':
                        print('Warning, have redundant %s in parallel with %s' % (name, Vname))

        return net, False

    def _simplify_series(self, cptnames=None, explain=False):

        net, changed = self._simplify_redundant_series(cptnames, explain)        
        net, changed2 = net._simplify_combine_series(cptnames, explain)
        return net, changed or changed2

    def _simplify_parallel(self, cptnames=None, explain=False):

        net, changed = self._simplify_redundant_parallel(cptnames, explain)
        net, changed2 = net._simplify_combine_parallel(cptnames, explain)
        return net, changed or changed2    

    def simplify_series(self, cptnames=None, explain=False, modify=True):

        net, changed = self._simplify_series(cptnames, explain)
        if not modify:
            return self
        return net

    def simplify_parallel(self, cptnames=None, explain=False, modify=True):

        net, changed = self._simplify_parallel(cptnames, explain)
        if not modify:
            return self        
        return net                

    def simplify(self, cptnames=None, passes=0, series=True,
                 parallel=True, explain=False, modify=True):

        # Perhaps use num cpts?
        if passes == 0:
            passes = 10                

        net = self
        for m in range(passes):

            net, series_changed = net._simplify_series(cptnames, explain)
            net, parallel_changed = net._simplify_parallel(cptnames, explain)
            if not series_changed and not parallel_changed:
                break
        if not modify:
            return self
            
        return net

    def check(self):
        """Check if network contains a loop of voltage sources or a cut set of current sources."""

        return self.simplify(explain=True, modify=False)
    
    def twoport(self, N1p, N1m, N2p, N2m, model='B'):
        """Create s-domain twoport model for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """

        from .twoport import TwoPortBModel, TwoPortHModel, TwoPortYModel, TwoPortZModel

        # TODO, generalise for not just s-domain.

        net = self.copy()
        if '0' not in net.nodes:
            net.add('W %s 0' % N1m)        

        if model == 'B':
            V2b = net.Voc(N2p, N2m)(s)
            I2b = net.Isc(N2p, N2m)(s)
            A = net.Aparams(N1p, N1m, N2p, N2m)
            return TwoPortBModel(A.Bparams, V2b, I2b)
        elif model == 'Z':
            V1 = net.Voc(N1p, N1m)(s)
            V2 = net.Voc(N2p, N2m)(s)
            Z = net.Zparams(N1p, N1m, N2p, N2m)            
            return TwoPortZModel(Z, V1, V2)
        elif model == 'Z':
            I1 = net.Isc(N1p, N1m)(s)
            I2 = net.Isc(N2p, N2m)(s)
            Z = net.Zparams(N1p, N1m, N2p, N2m)            
            return TwoPortYModel(Z.Y, I1, I2)        
        elif model == 'H':
            V1 = net.Voc(N1p, N1m)(s)
            I2 = net.Isc(N2p, N2m)(s)
            Z = net.Zparams(N1p, N1m, N2p, N2m)            
            return TwoPortHModel(Z.H, V1, I2)
        else:
            raise ValueError('Model %s unknown, must be B, H, Y, or Z' % model)

    @property
    def sch(self):
        """Generate schematic of subnetlist."""                

        if hasattr(self, '_sch'):
            return self._sch

        sch = Schematic(allow_anon=self.allow_anon)

        netlist = self.netlist()
        for net in netlist.split('\n'):
            sch.add(net)

        self._sch = sch
        return sch

    @property
    def sim(self):
        """Generate simulation object."""        

        if hasattr(self, '_sim'):
            return self._sim

        self._sim = Simulator(self)
        return self._sim
    
    @property
    def ss(self):
        """Generate state-space representation."""        

        if hasattr(self, '_ss'):
            return self._ss

        self._ss = StateSpace(self)
        return self._ss

    def replace(self, oldname, newname):
        """Replace component.
        
        For example,
        b = a.replace('C', 'W')
        c = a.replace('C1', 'C1 1 2')
        """        

        new = self._new()

        newparts = newname.split(' ')

        for cpt in self._elements.values():
            net = cpt._copy()            
            if cpt.name == oldname:
                if len(newparts) == 1:
                    # Just replace name of component
                    parts = net.split(' ')
                    parts[0] = newname
                    net = ' '.join(parts)
                else:
                    # Replace with new net
                    net = newname
            new._add(net)
        return new

    def subs(self, subs_dict):
        """Substitute values using dictionary of substitutions.
        
        For example, b = a.subs({'R1': 1e3, 'R2': 9e3})"""

        new = self._new()

        for cpt in self._elements.values():
            net = cpt._subs(subs_dict)
            new._add(net)
        return new                        
    
    def initialize(self, cct, time):
        """Set the initial values for this netlist based on the values
        computed for netlist cct at specified time."""        

        new = self._new()

        for cpt in self._elements.values():
            ic = 0
            if cpt.name in cct.reactances:
                if cpt.type == 'C':
                    ic = cct[cpt.name].v.remove_condition()(time)
                else:
                    ic = cct[cpt.name].i.remove_condition()(time) 
                    
            net = cpt._initialize(ic)
            new._add(net)
        return new                
    
    def pre_initial_model(self):
        """Generate circuit model for determining the pre-initial
        conditions."""

        new = self._new()

        for cpt in self._elements.values():
            net = cpt._pre_initial_model()
            new._add(net)
        return new        

    def r_model(self):
        """"Create resistive equivalent model using companion circuits.
        This is experimental!"""

        new = self._new()

        for cpt in self._elements.values():
            net = cpt._r_model()
            new._add(net)
        return new
    
    def s_model(self, var=s):
        """"Create s-domain model."""

        new = self._new()

        for cpt in self._elements.values():
            net = cpt._s_model(var)
            new._add(net)
        return new

    def ss_model(self):
        """"Create state-space model by replacing inductors
        with current sources and capacitors with voltage sources."""

        new = self._new()

        for cpt in self._elements.values():
            net = cpt._ss_model()
            new._add(net)
        return new

    def state_space_model(self):
        """"Create state-space model by replacing inductors
        with current sources and capacitors with voltage sources."""
        
        return self.ss_model()

    def ac_model(self, var=omega):
        """"Create AC model for specified angular frequency (default
        omega)."""
        
        return self.s_model(j * var)

    def noise_model(self, T='T'):
        """"Create noise model where resistors are converted into a series
        combination of an ideal resistor and a noise voltage
        source."""
        
        return self.noisy(T=T)
    
    def draw(self, filename=None, **kwargs):
        """Draw schematic of netlist.

        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        Note, if using Jupyter, then need to first issue command %matplotlib inline

        kwargs include:
           label_ids: True to show component ids
           label_values: True to display component values
           draw_nodes: True to show all nodes, False to show no nodes, 
             'primary' to show primary nodes,
             'connections' to show nodes that connect more than two components,
             'all' to show all nodes
           label_nodes: True to label all nodes, False to label no nodes, 
             'primary' to label primary nodes,
             'alpha' to label nodes starting with a letter,
             'pins' to label nodes that are pins on a chip,
             'all' to label all nodes
           style: 'american', 'british', or 'european'
           scale: schematic scale factor, default 1.0
           node_spacing: spacing between component nodes, default 2.0
           cpt_size: size of a component, default 1.5
           dpi: dots per inch for png files
           help_lines: distance between lines in grid, default 0.0 (disabled)
           debug: True to display debug information
        """

        cct = self
        if kwargs.pop('s_model', False):
            cct = cct.s_model()

        return cct.sch.draw(filename=filename, **kwargs)

    @property
    def is_causal(self):
        """Return True if all independent sources are causal and not an
        initial value problem (unless all the initial values are zero)."""
        return self.analysis['causal']        

    @property
    def is_dc(self):
        """Return True if all independent sources are DC and not an
        initial value problem.  The initial value problem may collapse
        to a DC problem but we cannot prove this yet."""
        return self.analysis['dc']

    @property
    def is_ac(self):
        """Return True if all independent sources are AC and not an
        initial value problem."""
        return self.analysis['ac']

    @property
    def is_superposition(self):
        """Return True if netlist needs to be solved using multiple approaches,
        e.g., ac and dc"""
        
        return ((self.has_ac and self.has_dc) or
                (self.has_ac and self.has_transient) or
                (self.has_dc and self.has_transient))

    @property
    def has_dc(self):
        """Return True if any independent source has a DC component."""
        return self.analysis['has_dc']

    @property
    def has_ac(self):
        """Return True if any independent source has an AC component."""
        return self.analysis['has_ac']

    @property
    def has_s_transient(self):
        """Return True if any independent source has a transient component defined in s-domain."""
        return self.analysis['has_s']    

    @property    
    def has_transient(self):
        """Return True if any independent source has a transient component."""
        return self.analysis['has_transient']    
    
    @property
    def zeroic(self):
        """Return True if the initial conditions for all components are zero."""
        return self.analysis['zeroic']

    @property
    def is_ivp(self):
        """Return True for an initial value problem.  This is True if any
        component that allows initial conditions has them explicitly
        defined.

        """
        return self.analysis['ivp']

    @property
    def is_time_domain(self):
        """Return True if can analyse in time domain."""
        return self.analysis['time_domain']

    @property
    def missing_ic(self):
        """Return components that allow initial conditions but do not have
        them explicitly defined."""

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.has_ic is False)

    @property
    def sources(self):
        """Return dictionary of all sources (this does not include
        implicit sources due to initial conditions)."""

        return self.dependent_sources + self.independent_sources

    @property
    def reactances(self):
        """Return dictionary of reactances."""

        return self.analysis['reactances']

    @property
    def ics(self):
        """Return dictionary of components with initial conditions."""

        return self.analysis['ics']    
    
    @property
    def independent_sources(self):
        """Return dictionary of independent sources (this does not include
        implicit sources due to initial conditions)."""

        return self.analysis['independent_sources']

    @property
    def dependent_sources(self):
        """Return dictionary of dependent sources."""

        return self.analysis['dependent_sources']            

    def independent_source_groups(self, transform=False):
        """Return dictionary of source groups.  Each group is a list of
        sourcenames that can be analysed at the same time.  Noise
        sources have separate groups since they are assumed to be
        uncorrelated.  The source groups are keyed by the names of
        superposition decomposition keys: 'dc', 's', 't', 'n*', and omega,
        where omega is an expression for the angular frequency of a phasor,
        and `n*` is a noise identifier.
        
        If transform is False, the returned keys are 'dc', 't', 's', and 'n*'.

        If transform is True, the returned keys are 'dc', 's', omega,
        and 'n*'.  Note, there is no time-domain component.  Note,
        after transformation, a source can appear in multiple groups.
        For example, if a voltage source V1 has a value 10 + 5 *
        cos(omega * t), V1 will be added to the dc and omega groups.

        """

        groups = {}
        for eltname, elt in self.elements.items():
            if not elt.independent_source:
                continue
            cpt = elt.cpt
            if cpt.is_voltage_source:
                Voc = cpt.Voc
                if transform:
                    Voc = Voc.decompose()
                cpt_kinds = Voc.keys()
            else:
                Isc = cpt.Isc
                if transform:
                    Isc = Isc.decompose()                
                cpt_kinds = Isc.keys()
                
            for cpt_kind in cpt_kinds:
                if cpt_kind not in groups:
                    groups[cpt_kind] = []
                groups[cpt_kind].append(eltname)

        return groups    

    @property
    def control_sources(self):
        """Return list of voltage sources required to specify control
        current for CCVS and CCCS components."""
        return self.analysis['control_sources']        

    @property
    def analysis(self):

        if not hasattr(self, '_analysis'):
            self._analysis = self.analyse()

        return self._analysis

    def analyse(self):

        has_ic = False
        zeroic = True
        has_s = False
        has_transient = False        
        ac_count = 0
        dc_count = 0
        causal = True
        reactive = False
        independent_sources = []
        dependent_sources = []        
        control_sources = []
        reactances = []
        ics = []
        
        for eltname, elt in self.elements.items():
            if elt.need_control_current:
                control_sources.append(elt.args[0])
            if elt.has_ic is not None:
                if elt.has_ic:
                    has_ic = True
                    ics.append(eltname)
                if not elt.zeroic:
                    zeroic = False
            if elt.independent_source:
                independent_sources.append(eltname)
                if elt.has_s_transient:
                    has_s = True
                if elt.has_transient:
                    has_transient = True                    
                if elt.is_ac:
                    ac_count += 1
                if elt.is_dc:
                    dc_count += 1
                if not elt.is_causal:
                    causal = False
            if elt.dependent_source:
                dependent_sources.append(eltname)
            if elt.reactive:
                reactive = True
                reactances.append(eltname)                

        num_sources = len(independent_sources)
                    
        analysis = {} 
        analysis['zeroic'] = zeroic
        analysis['has_ic'] = has_ic        
        analysis['ivp'] = has_ic
        analysis['has_dc'] = dc_count > 0
        analysis['has_ac'] = ac_count > 0        
        analysis['has_s'] = has_s
        analysis['has_transient'] = has_transient
        analysis['reactances'] = reactances
        analysis['ics'] = ics         
        analysis['dependent_sources'] = dependent_sources        
        analysis['independent_sources'] = independent_sources
        analysis['control_sources'] = control_sources        
        analysis['ac'] = ac_count > 0 and (num_sources == ac_count) and not has_ic
        analysis['dc'] = dc_count > 0 and (num_sources == dc_count) and not has_ic
        analysis['causal'] = causal and zeroic
        analysis['time_domain'] = not reactive and not has_s

        if not reactive and has_ic:
            raise ValueError('Non-reactive component with initial conditions')
        return analysis

    def describe(self):
        """Print a message describing how circuit is solved."""

        def describe_sources(sources):
            sources_string = ', '.join(sources)
            if len(sources) == 1:
                return 'source %s' % sources_string
            return 'sources %s' % sources_string

        def describe_analysis(method, sources):
            return '%s analysis is used for %s.' % (method,
                                                    describe_sources(sources))

        groups = self.independent_source_groups(transform=not
                                                self.is_time_domain)

        if groups == {}:
            print('There are no non-zero independent sources so everything is zero.')
            return

        if self.is_ivp:
            print('This has initial conditions for %s so is an initial value '
                  'problem solved in the s-domain using Laplace transforms.' %
                  ', '.join(self.ics))
            return

        if len(groups) > 1:
            print('This is solved using superposition.')
        for kind, sources in groups.items():
            if not isinstance(kind, str):
                print(describe_analysis('Phasor', sources))
            elif kind[0] == 'n':
                print(describe_analysis('Noise', sources))
            elif kind == 'dc':
                print(describe_analysis('DC', sources))
            elif kind == 's':
                print(describe_analysis('Laplace', sources))
            elif kind in ('t', 'time'):
                print(describe_analysis('Time-domain', sources))

    def Vname(self, name):
        return Vname(name, self.kind)

    def Iname(self, name):
        return Iname(name, self.kind)    


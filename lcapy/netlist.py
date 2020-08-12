"""This module provides support for the common aspects of Circuit and
Network classes.

Copyright 2014--2020 Michael Hayes, UCECE

"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from .sexpr import Hs, Zs, Ys
from .symbols import j, s, omega
from .voltage import Voltage, Vname
from .current import Current, Iname
from .schematic import Schematic
from .mna import MNAMixin, Nodedict, Branchdict
from .statespace import StateSpace
from .netfile import NetfileMixin
from .expr import Expr
from .state import state
from .attrdict import AttrDict
from .immitance import ImmitanceMixin
from .impedance import Impedance
from .admittance import Admittance
from . import mnacpts
from copy import copy
from collections import OrderedDict


class Node(ImmitanceMixin):

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

    @property
    def admittance(self):
        """Driving-point admittance between node and ground."""

        return self.cct.admittance(self.name, '0')

    @property
    def impedance(self):
        """Driving-point impedance between node and ground."""        

        return self.cct.impedance(self.name, '0')        

    def append(self, cpt):

        if cpt.type in ('P', ):
            self.port = True

        self.list.append(cpt)

    def oneport(self, node=0):
        """Create oneport object with respect to specified node
        (default ground)."""

        return self.cct.oneport(self.name, node)        
        
    def thevenin(self, node=0):
        """Create Thevenin oneport object with respect to specified
        node (default ground)."""
        
        return self.cct.thevenin(self.name, node)

    def norton(self, node=0):
        """Create Norton oneport object with respect to specified node
        (default ground)."""
        
        return self.cct.norton(self.name, node)    
        

class NetlistNamespace(object):
    """This class allows elements, nodes, or other namespaces
    in a heirachical namespace to be accessed by name, via __getattr__.
    
    For example,

    >>> a = Circuit('''
    ... b.L1 1 2'''
    >>> a.b.L1.v

    """

    def __init__(self, namespace, netlist):

        # The namespace name, e.g., a or a.b
        self.namespace = namespace
        self._netlist = netlist
        # The children namespaces
        self.namespaces = {}
        self.allow_anon = False

    def __getitem__(self, name):
        """Return element or node by name."""

        netlist = self._netlist
        
        # If name is an integer, convert to a string.
        if isinstance(name, int):
            name = '%d' % name

        name = self.namespace + '.' + name
            
        if name in netlist.nodes:
            return netlist.nodes[name]

        if name in netlist._elements:
            return netlist._elements[name]

        # Try first anonymous name.
        if name + 'anon1' in netlist._elements:
            return netlist._elements[name + 'anon1']

        if name in self.namespaces:
            return self.namespaces[name]
        
        raise AttributeError('Unknown element or node name %s' % name)

    def __getattr__(self, attr):
        """Return element, node, or another NetlistNamespace object by name.
        This gets called if there is no explicit attribute attr for
        this instance.  This is primarily for accessing elements and
        non-numerical node names.  It also gets called if the called
        attr throws an AttributeError exception.  The annoying thing
        is that hasattr uses getattr and checks for an exception.

        """

        return self.__getitem__(attr)

    def netlist(self):
        """Return the current netlist for this namespace."""

        nlist = self._netlist
        
        return '\n'.join([str(cpt) for cpt in nlist._elements.values() if str(cpt).startswith(self.namespace)])

    def __repr__(self):
        
        return self.netlist()

    @property
    def sch(self):
        """Generate schematic for this namespace."""        

        if hasattr(self, '_sch'):
            return self._sch

        sch = Schematic(allow_anon=self.allow_anon)

        netlist = self._netlist.netlist()
        for net in netlist.split('\n'):
            if net.startswith(self.namespace):
                sch.add(net)

        self._sch = sch
        return sch

    def draw(self, filename=None, **kwargs):
        """Draw schematic of subnetlist.

        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.
        
        Note, if using Jupyter, then need to first issue command %matplotlib inline

        kwargs include:
           label_ids: True to show component ids
           label_values: True to display component values
           draw_nodes:
              True to show all nodes,
              False or 'none' to show no nodes, 
              'primary' to show primary nodes,
              'connections' to show nodes that connect more than two components,
              'all' to show all nodes
           label_nodes:
               True to label all nodes,
               False or 'none' to label no nodes, 
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
            cct = cct._s_model()

        return cct.sch.draw(filename=filename, **kwargs)


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

        for node in cpt.nodes:
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
        return self        

    @property
    def equipotential_nodes(self):
        """Determine nodes connected by wires that are of the same potential.
        This returns a dictionary keyed by the unique node names with
        values being lists of nodes of the same potential."""

        enodes = {}
        for key in self.nodes.keys():
            enodes[key] = [key]

        # Then augment with nodes connected by wires.
        for m, elt in enumerate(self.elements.values()):
            if elt.type not in ('W', ):
                continue

            n1, n2 = elt.nodes

            for key1, nodes in enodes.items():
                if n1 in nodes:
                    break

            for key2, nodes in enodes.items():
                if n2 in nodes:
                    break

            if key1 != key2:
                enodes[key1].extend(enodes.pop(key2))

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
    
    def augment_node_map(self, node_map={}):

        """Create a mapping dict for all nodes."""

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
    
    def renumber(self, node_map={}):
        """Renumber nodes using specified node_map.  If node_map not specified
        then a mapping is created."""

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
                Np, Nm = cpt.nodes[0:2]
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

        return Admittance(If.laplace(), kind=self.kind)

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

        return Impedance(Vf.laplace(), kind=self.kind)

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
        
    def transfer(self, N1p, N1m, N2p, N2m):
        """Create s-domain voltage transfer function V2(s) / V1(s) where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed and initial conditions
        are ignored.  Since the result is causal, the frequency response
        can be found by substituting j * omega for s."""

        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)
        
        new = self.kill()
        if '0' not in new.nodes:
            new.add('W %s 0' % N1m)
        
        new._add('V1_ %s %s {DiracDelta(t)}' % (N1p, N1m))

        V2 = new.Voc(N2p, N2m)
        V1 = new.V1_.V

        return Hs(V2.laplace() / V1.laplace(), causal=True)

    def Amatrix(self, N1p, N1m, N2p, N2m):
        """Create A matrix from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
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
            A11 = Hs(net.V1_.V(s) / net.Voc(N2p, N2m)(s))

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = Zs(net.V1_.V(s) / net.Isc(N2p, N2m)(s))

            net.remove('V1_')

            net.add('I1_ %s %s {DiracDelta(t)}' % (N1p, N1m))

            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure I2 with port 2 open-circuit
            try:
                A21 = Ys(-1 / net.Voc(N2p, N2m)(s))                
            except ValueError:
                # It is likely there is an open-circuit.                
                net2 = net.copy()
                net2.add('W %s %s' % (N2p, N2m))
                A21 = Ys(-net2.I1_.I(s) / net2.Voc(N2p, N2m)(s))
                A21 = 0                

            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure I2 with port 2 short-circuit
            A22 = Hs(1 / net.Isc(N2p, N2m)(s))

            net.remove('I1_')
            A = AMatrix(((A11, A12), (A21, A22)))
            return A

        except ValueError:
            raise ValueError('Cannot create A matrix')


    def Zmatrix(self, N1p, N1m, N2p, N2m):
        """Create Z matrix from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """

        from .twoport import ZMatrix

        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)        
        net = self.kill()
        if '0' not in net.nodes:
            net.add('W %s 0' % N1m)

        try:
            net.add('I1_ %s %s {DiracDelta(t)}' % (N1p, N1m))

            # Z11 = V1 / I1 with I2 = 0
            # Apply I1 and measure V1 with port 2 open-circuit
            Z11 = Zs(net.Voc(N1p, N1m)(s))

            # Z21 = V2 / I1 with I2 = 0
            # Apply I1 and measure V2 with port 2 open-circuit
            Z21 = Zs(net.Voc(N2p, N2m)(s))          

            net.remove('I1_')

            net.add('I2_ %s %s {DiracDelta(t)}' % (N2p, N2m))

            # Z12 = V1 / I2 with I1 = 0
            # Apply I2 and measure V1 with port 1 open-circuit
            Z12 = Zs(net.Voc(N1p, N1m)(s))

            # Z22 = V2 / I2 with I1 = 0
            # Apply I2 and measure V2 with port 1 open-circuit
            Z22 = Zs(net.Voc(N2p, N2m)(s))          

            net.remove('I2_')            

            Z = ZMatrix(((Z11, Z12), (Z21, Z22)))
            return Z

        except ValueError:
            raise ValueError('Cannot create Z matrix')        

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

    def _noisy(self, resistornames):

        new = self._new()

        for cpt in self._elements.values():
            if cpt.name in resistornames:
                net = cpt._noisy()
            else:
                net = cpt._copy()
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

    def twoport(self, N1p, N1m, N2p, N2m, model='B'):
        """Create s-domain twoport model from network, where:
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
            A = net.Amatrix(N1p, N1m, N2p, N2m)
            return TwoPortBModel(A.B, V2b, I2b)
        elif model == 'Z':
            V1 = net.Voc(N1p, N1m)(s)
            V2 = net.Voc(N2p, N2m)(s)
            Z = net.Zmatrix(N1p, N1m, N2p, N2m)            
            return TwoPortZModel(Z, V1, V2)
        elif model == 'Z':
            I1 = net.Isc(N1p, N1m)(s)
            I2 = net.Isc(N2p, N2m)(s)
            Z = net.Zmatrix(N1p, N1m, N2p, N2m)            
            return TwoPortYModel(Z.Y, I1, I2)        
        elif model == 'H':
            V1 = net.Voc(N1p, N1m)(s)
            I2 = net.Isc(N2p, N2m)(s)
            Z = net.Zmatrix(N1p, N1m, N2p, N2m)            
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
    def ss(self):
        """Generate state-space representation."""        

        if hasattr(self, '_ss'):
            return self._ss

        self._ss = StateSpace(self)
        return self._ss

    def subs(self, subs_dict):
        """Substitute values using dictionary of subsitutions."""        

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

    def noise_model(self):
        """"Create noise model where resistors are converted into a series
        combination of an ideal resistor and a noise voltage
        source."""
        
        return self.noisy()
    
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
    def has_dc(self):
        """Return True if any independent source has a DC component."""
        return self.analysis['has_dc']

    @property
    def has_ac(self):
        """Return True if any independent source has an AC component."""
        return self.analysis['has_ac']

    @property
    def has_s(self):
        """Return True if any independent source has an s-domain component."""
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

        return dict((key, cpt) for key, cpt in self.elements.items() if cpt.hasic is False)

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
        uncorrelated.

        If transform is True, the source values are decomposed into
        the different transform domains to determine which domains are
        required.  In this case, a source can appear in multiple
        groups.  For example, if a voltage source V1 has a value 10 +
        5 * cos(omega * t), V1 will be added to the dc and ac groups.

        """

        groups = {}
        for key, elt in self.elements.items():
            if not elt.independent_source:
                continue
            cpt = elt.cpt
            if cpt.voltage_source:
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
                groups[cpt_kind].append(key)

        return groups    

    @property
    def control_sources(self):
        """Return dictionary of voltage sources required to specify control
        current for CCVS and CCCS components."""
        return self.analysis['control_sources']        

    @property
    def analysis(self):

        if not hasattr(self, '_analysis'):
            self._analysis = self.analyse()

        return self._analysis

    def analyse(self):

        hasic = False
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
        
        for key, elt in self.elements.items():
            if elt.need_control_current:
                control_sources.append(elt.args[0])
            if elt.hasic is not None:
                if elt.hasic:
                    hasic = True
                if not elt.zeroic:
                    zeroic = False
            if elt.independent_source:
                independent_sources.append(key)
                if elt.has_s:
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
                dependent_sources.append(key)
            if elt.reactive:
                reactive = True
                reactances.append(key)                

        num_sources = len(independent_sources)
                    
        analysis = {} 
        analysis['zeroic'] = zeroic
        analysis['has_ic'] = hasic        
        analysis['ivp'] = hasic
        analysis['has_dc'] = dc_count > 0
        analysis['has_ac'] = ac_count > 0        
        analysis['has_s'] = has_s
        analysis['has_transient'] = has_transient
        analysis['reactances'] = reactances 
        analysis['dependent_sources'] = dependent_sources        
        analysis['independent_sources'] = independent_sources
        analysis['control_sources'] = control_sources        
        analysis['ac'] = ac_count > 0 and (num_sources == ac_count) and not hasic
        analysis['dc'] = dc_count > 0 and (num_sources == dc_count) and not hasic
        analysis['causal'] = causal and zeroic
        analysis['time_domain'] = not reactive and not has_s

        if not reactive and hasic:
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
            print('This has initial conditions so is an initial value problem '
                  'solved in the s-domain using Laplace transforms.')
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
            elif kind == 'time':
                print(describe_analysis('Time-domain', sources))

    def Vname(self, name):
        return Vname(name, self.kind)

    def Iname(self, name):
        return Iname(name, self.kind)    

                
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
        if isinstance(key, Expr):
            key = key.expr
        return super(Transformdomains, self).__getitem__(key)

    
class Netlist(NetlistMixin, NetfileMixin):
    """This class handles a generic netlist with multiple sources.
    During analysis, subnetlists are created for each source kind (dc,
    ac, transient, etc).  Since linearity is assumed, superposition is
    employed.

    """

    def __init__(self, filename=None, context=None, allow_anon=False):

        super (Netlist, self).__init__(filename, context, allow_anon=allow_anon)
        self._invalidate()
        self.kind = 'super'

    def _invalidate(self):

        for attr in ('_sch', '_sub', '_Vdict', '_Idict', '_analysis',
                     '_node_map', '_ss', '_node_list', '_branch_list'):
            try:
                delattr(self, attr)
            except:
                pass

    def _sub_make(self):

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

        elif self.is_time_domain:
            groups = self.independent_source_groups()            
            newgroups = {'time' : []}
            for key, sources in groups.items():
                if isinstance(key, str) and key[0] == 'n':
                    newgroups[key] = sources
                else:
                    newgroups['time'] += sources
            groups = newgroups

        else:
            groups = self.independent_source_groups(transform=True)        
        
        self._sub = Transformdomains()

        for kind, sources in groups.items():
            self._sub[kind] = SubNetlist(self, kind)

        return self._sub
        
    @property
    def sub(self):
        """Return dictionary of subnetlists keyed by transform domain kind.
        Note, the subnetlists are not created until a specific one is
        selected.

        """

        if hasattr(self, '_sub'):
            return self._sub

        return self._sub_make()

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
                    result[node] = Voltage()
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
                    result[node] = Current()
                result[node].add(value)
        self._Idict = result                    
        return result    

    def get_I(self, name):
        """Current through component"""

        result = Current()
        for sub in self.sub.values():
            I = sub.get_I(name)
            result.add(I)
        result = result.canonical()            
        return result

    def get_i(self, name):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def _get_Vd(self, Np, Nm=None):
        """This does not check nodes."""
        
        result = Voltage()
        for sub in self.sub.values():
            Vd = sub.get_Vd(Np, Nm)
            result.add(Vd)
        result = result.canonical()
        return result

    def get_Vd(self, Np, Nm=None):
        """Voltage drop between nodes"""

        Np, Nm = self._check_nodes(Np, Nm)
        return self._get_Vd(Np, Nm)

    def get_vd(self, Np, Nm=None):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(Np, Nm).time()

    def dc(self):
        """Return subnetlist for dc components of independent sources.

        See also, ac, transient, laplace.
        """
        return SubNetlist(self, 'dc')

    def ac(self):
        """Return subnetlist for ac components of independent sources
        for angular frequency omega.

        See also, dc, transient, laplace.
        """
        # Could look at all the ac frequencies and if there is only
        # one use that?  If have multiple ac frequencies should issue
        # warning.
        return SubNetlist(self, omega)    

    def transient(self):
        """Return subnetlist for transient components of independent
        sources.  Note, unlike the similar laplace method, dc and ac 
        components are ignored.

        See also, dc, ac, laplace.

        """        
        return SubNetlist(self, 's')

    def laplace(self):
        """Return subnetlist for Laplace representations of independent
        source values.

        See also, dc, ac, transient.
        
        """        
        return SubNetlist(self, 'laplace')    
    
    
class SubNetlist(NetlistMixin, MNAMixin):
    """This is a representation of a netlist for a particular
    transformation domain, such as ac, dc, transient, or noise."""

    def __new__(cls, netlist, kind):

        obj = netlist.select(kind=kind)
        # Need own context to avoid conflicts with Vn1 and Vn1(s), etc.
        obj.context = state.new_context()
        obj.kind = kind
        obj.__class__ = cls
        obj._analysis = obj.analyse()
        return obj

    def __init__(cls, netlist, kind):
        """ kind can be 't', 'dc', 'ac', 's', 'time', 'ivp', 'n*', omega, 
        or an integer"""
        
        if kind == omega:
            return
        if not isinstance(kind, str):
            return
        if kind[0] == 'n':
            return
        kinds = ('t', 'dc', 'ac', 's', 'time', 'ivp', 'laplace')
        if kind not in kinds:
            raise ValueError('Expected one of %s for kind, got %s' %
                             (', '.join(kinds), kind))

    def get_I(self, name):
        """Current through component"""

        self._solve()
        return self._Idict[name].canonical()

    def get_i(self, name):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def get_Vd(self, Np, Nm=None):
        """Voltage drop between nodes"""

        self._solve()
        return (self._Vdict[Np] - self._Vdict[Nm]).canonical()

    def get_vd(self, Np, Nm=None):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(Np, Nm).time()

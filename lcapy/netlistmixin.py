"""This module provides the NetlistMixin class.  This is used for
Netlist and SubNetlist.

Copyright 2020--2023 Michael Hayes, UCECE

"""

from .expr import expr
from .analysis import Analysis
from .components import Components
from .mnacpts import Cpt
from .impedance import impedance
from .admittance import admittance
from .equipotentialnodes import EquipotentialNodes
from .node import Node
from .state import state
from .symbols import j, s, omega
from .attrdict import AttrDict
from .netfile import NetfileMixin
from .statespace import StateSpace
from .voltage import Vname
from .current import Iname, current
from .transfer import transfer
from .simulator import Simulator
from .netlistnamespace import NetlistNamespace
from .matrix import Matrix
from .utils import isiterable
from .texpr import t
from .deprecation import LcapyDeprecationWarning

from . import mnacpts
from collections import OrderedDict
from warnings import warn


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

        import pdb
        pdb.set_trace()
        return self

    @property
    def cpts(self):
        """Return list of component names."""

        # Perhaps should prune wires, open-circuits, etc. ?
        return list(self._elements.keys())

    def _dummy_node(self):
        """Create a dummy node name."""

        return '_' + self._make_anon_name('node')

    def _add_ground(self, node):

        if '0' not in self.nodes:
            self.add('W %s 0' % node)

    def _add_test_voltage_source(self, Np, Nm):

        self._add('V? %s %s {DiracDelta(t)}' % (Np, Nm))
        return self.last_added()

    def _add_test_current_source(self, Np, Nm):

        self._add('I? %s %s {DiracDelta(t)}' % (Np, Nm))
        return self.last_added()

    def apply_test_voltage_source(self, Np, Nm=None):
        """This copies the netlist, kills all the sources, and applies a Dirac
        delta test voltage source across the specified nodes.  If the
        netlist is not connected to ground, the negative specified
        node is connected to ground.  The new netlist is returned."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.kill()
        new._add_ground(Nm)
        new._add_test_voltage_source(Np, Nm)
        return new

    def apply_test_current_source(self, Np, Nm=None):
        """This copies the netlist, kills all the sources, and applies a Dirac
        delta test current source across the specified nodes.  If the
        netlist is not connected to ground, the negative specified
        node is connected to ground.  The new netlist is returned."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.kill()
        new._add_ground(Nm)
        new._add_test_current_source(Np, Nm)
        return new

    @property
    def params(self):
        """Return list of symbols used as arguments in the netlist."""

        symbols = self.symbols
        params = []
        for elt in self.elements.values():
            for arg in elt.args:
                if arg in symbols and arg not in params:
                    params.append(arg)
        return params

    @property
    def symbols(self):
        """Return dictionary of symbols defined in the netlist."""

        return self.context.symbols

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
            warn('Overriding component %s' % cpt.name)
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

    def expand(self):
        """Expand the netlist, replacing complicated components with simpler
        components."""

        new = self._new()

        for cpt in self._elements.values():
            new._add(cpt._expand())
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
        values being lists of nodes of the same potential.

        Primary nodes (without underscore) are chosen for the key
        if possible."""

        enodes = EquipotentialNodes()
        enodes.add(self.nodes.keys())

        # Then augment with nodes connected by wires.
        for m, elt in enumerate(self.elements.values()):
            if elt.nosim:
                continue
            if elt.type == 'W':
                enodes.add_wire(*elt.nodenames)
            elif elt.type.startswith('TL'):
                enodes.add_wire(elt.nodenames[1], elt.nodenames[3])
            elif elt.type.startswith('TP'):
                enodes.add_wire(elt.nodenames[1], elt.nodenames[3])
                warn("Assuming V2' = V1' for %s" % elt.name)
            else:
                for connections in elt.equipotential_nodes:
                    enodes.add_wires([elt.name + '.' + n for n in connections])

        # Alter keys to avoid underscore and to ensure that have a '0'
        # key if possible.
        enodes2 = {}
        for key, nodes in enodes.items():
            # Sort so that lowest primary nodes come first
            nodes = sorted(nodes, key=lambda item: ('_' in item, item))
            if '0' in nodes:
                newkey = '0'
            else:
                newkey = nodes[0]
            enodes2[newkey] = nodes

        return enodes2

    def unconnected_nodes(self):
        """Return list of node names that are not connected."""

        return [node.name for node in self.nodes.values() if node.count == 1]

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

    def annotate_currents(self, cpts=None, domainvar=None, flow=False,
                          eng_format=True, evalf=True, num_digits=3,
                          show_units=True, pos=''):
        """Annotate specified list of component names `cpts` with current (or
        flow).

        `domainvar` specifies the domain to calculate the voltages for
        (e.g., `t` for time-domain, `s` for Laplace-domain)

        `flow` (default False) if True annotates current as a flow

        `eng_format` (default True) if True use engineering format if
        the current is a number, e.g., 100\,mV instead of 0.1\,V

        `evalf` (default True) if True prints floating point
        numbers as decimals otherwise they are shown as rationals

        `show_units` (default True) if True applies the units (e.g.,
        V for volts)

        `pos` specifies where to position the labels (see docs)
        """

        if cpts is None:
            cpts = []
            for elt in self._elements.values():
                if (elt.is_resistor or elt.is_capacitor or
                        elt.is_inductor or elt.is_voltage_source):
                    cpts.append(elt.name)

        label = ('f' if flow else 'i') + pos

        if domainvar is None:
            domainvar = t

        new = self._new()
        for cpt in self._elements.values():
            net = cpt._copy()
            if cpt.name in cpts:
                I = cpt.I(domainvar)
                if evalf:
                    I = I.evalf(num_digits)
                net += ', ' if ';' in net else '; '
                net += ', %s={$%s$}' % (label, I.latex_with_units(
                    eng_format=eng_format, evalf=evalf, num_digits=num_digits, show_units=show_units))
            new.add(net)
        return new

    def annotate_voltages(self, cpts, domainvar=None,
                          eng_format=True, evalf=True, num_digits=3,
                          show_units=True, pos=''):
        """Annotate specified list of component names `cpts` with voltage.

        `domainvar` specifies the domain to calculate the voltages for
        (e.g., `t` for time-domain, `s` for Laplace-domain)

        `pos` specifies where to position the labels, see docs

        `eng_format` (default True) if True use engineering format if
        the voltage is a number, e.g., 100\,mV instead of 0.1\,V

        `evalf` (default True) if True prints floating point
        numbers as decimals otherwise they are shown as rationals

        `show_units` (default True) if True applies the units (e.g.,
        V for volts)

        `pos` specifies where to position the labels (see docs)
        """

        if cpts is None:
            cpts = []
            for elt in self._elements.values():
                if (elt.is_resistor or elt.is_capacitor or
                        elt.is_inductor or elt.is_current_source):
                    cpts.append(elt.name)

        if domainvar is None:
            domainvar = t

        new = self._new()
        for cpt in self._elements.values():
            net = cpt._copy()
            if cpt.name in cpts:
                V = cpt.V(domainvar)
                if evalf:
                    V = V.evalf(num_digits)
                net += ', ' if ';' in net else '; '
                net += 'v%s={$%s$}' % (pos, V.latex_with_units(eng_format=eng_format,
                                       evalf=evalf, num_digits=num_digits, show_units=show_units))
            new.add(net)
        return new

    def annotate_node_voltages(self, nodes=None, domainvar=None,
                               label_voltages=False, eng_format=True,
                               evalf=True, num_digits=3,
                               show_units=True, anchor='south west'):
        """Create a new netlist with the node voltages annotated.  This is
        useful for drawing a schematic with the node voltages shown.
        For example,

        `cct.annotate_node_voltages((1, 2, 3)).draw()`

        `nodes` is a list of the nodes to annotate or `None` for all.

        `domainvar` specifies the domain to calculate the voltages for
        (e.g., `t` for time-domain, `s` for Laplace-domain)

        `label_voltages` (default False) if True prefixes the
        annotation with V1= for node 1, etc.

        `eng_format` (default True) if True use engineering format if
        the voltage is a number, e.g., 100\,mV instead of 0.1\,V

        `evalf` (default True) if True prints floating point
        numbers as decimals otherwise they are shown as rationals

        `num_digits` (default 3) specfies the number of digits to print
        for floating point numbers

        `show_units` (default True) if True applies the units (e.g.,
        V for volts)

        `anchor` (default 'south west') specifies the position of the
        voltage label

        """

        if nodes is None:
            nodes = self.node_list
        elif not isiterable(nodes):
            nodes = (nodes, )

        if domainvar is None:
            domainvar = t

        new = self.copy()
        for node in nodes:
            v = self[node].V(domainvar)
            if evalf:
                v = v.evalf(num_digits)

            vstr = '%s' % v.latex_with_units(eng_format=eng_format, evalf=evalf,
                                             num_digits=num_digits,
                                             show_units=show_units)

            if label_voltages:
                vstr = 'V_{%s}=' % node + vstr

            new.add('A%s %s; l={%s}, anchor=%s' % (node, node, vstr, anchor))
        return new

    def annotate_voltage(self, cpts, domainvar=None, pos=''):
        LcapyDeprecationWarning(
            feature="annotate_voltage",
            useinstead="annotate_voltages",
            deprecated_since_version="1.0",
            issue=61
        ).warn()
        return self.annotate_voltages(cpts, domainvar=domainvar,
                                      pos=pos, evalf=True, num_digits=3)

    def annotate_current(self, cpts, domainvar=None, flow=False, pos=''):
        LcapyDeprecationWarning(
            feature="annotate_current",
            useinstead="annotate_currents",
            deprecated_since_version="1.0",
            issue=61
        ).warn()
        return self.annotate_currents(cpts, domainvar=domainvar, flow=flow,
                                      pos=pos, evalf=True, num_digits=3)

    def augment_node_map(self, node_map=None):
        """Create a mapping dict for all nodes."""

        if node_map is None:
            node_map = {}

        # Add nodes with names like U1.vdd so they don't get renamed.
        ignore_nodes = []
        for node in self.nodes.keys():
            if node in node_map:
                continue
            parts = node.split('.')
            if len(parts) < 2:
                continue
            if parts[-2] in self._elements:
                node_map[node] = node
                ignore_nodes.append(node)

        # It would be desirable to renumber the nodes say from left to
        # right and top to bottom.  The schematic drawing algorithms
        # could help with this since they figure out the node
        # placement.

        enodes = self.equipotential_nodes

        # Rewrite the equipotential nodes removing nodes that need to
        # be ignored.
        enodes2 = {}
        for key, nodes in enodes.items():
            newnodes = []
            for node in nodes:
                if node not in ignore_nodes:
                    newnodes.append(node)

            if key in ignore_nodes:
                key = newnodes[0]
            enodes2[key] = newnodes
        enodes = enodes2

        if '0' in self.nodes and '0' not in node_map:
            node_map['0'] = '0'

        numbers = []
        for m in range(len(enodes)):
            numbers.append('%s' % (m + 1))

        # Check if user has supplied an unknown node in node_map.
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
                    raise ValueError(
                        'Cannot rename two nodes of same potential')
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
        ['C1', 'V1', 'R1', 'R2']."""

        if hasattr(self, '_branch_list'):
            return self._branch_list

        self._branch_list = []
        for key, elt in self.elements.items():
            if elt.type not in ('W', 'O', 'P', 'K', 'XX'):
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

    def _parse_node_args2(self, Np, Nm=None):

        if Nm is None:
            cpt = self[Np]
            if isinstance(cpt, Node):
                Np, Nm = cpt.name, 0
            else:
                Np, Nm = cpt.nodenames[0:2]
        return Np, Nm

    def _parse_node_args4(self, N1p, N1m, N2p, N2m, name):

        if N2p is None and N2m is None:

            arg1, arg2 = N1p, N1m

            if isinstance(arg1, tuple):
                N1p, N1m = arg1
                # TODO: check if there is a voltage source across these nodes
                # that will short out the applied source.
            else:
                try:
                    arg1 = self.elements[arg1]
                except:
                    pass
                if arg1 not in self.elements.values():
                    raise ValueError('Unknown component %s' % arg1)
                if arg1.is_voltage_source:
                    # The killed voltage source will short the applied signal.
                    raise ValueError(
                        "Cannot determine transfer function across voltage source %s; you will need to remove it, e.g., new = cct.remove('%s')" % (arg1, arg1))
                N1p, N1m = [n.name for n in arg1.nodes[0:2]]

            if isinstance(arg2, tuple):
                N2p, N2m = arg2
            else:
                try:
                    arg2 = self.elements[arg2]
                except:
                    pass

                if arg2 not in self.elements.values():
                    raise ValueError('Unknown component %s' % arg2)
                N2p, N2m = [n.name for n in arg2.nodes[0:2]]

        elif N2p is None or N2m is None:
            raise ValueError('Expecting %s(cpt1, cpt2), %s(cpt1, (N2p, N2m), %s((N1p, N1m), cpt2), or %s(N1p, N1m, N2p, N2m)' % (
                name, name, name, name))

        return N1p, N1m, N2p, N2m

    def Voc(self, Np, Nm=None, **kwargs):
        """Return open-circuit transform-domain voltage between nodes Np and
        Nm."""

        return self.get_Vd(Np, Nm, **kwargs)

    def voc(self, Np, Nm=None):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).time()

    def Isc(self, Np, Nm=None, **kwargs):
        """Return short-circuit transform-domain current between nodes Np and
        Nm."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.copy()
        if new.is_causal:
            new.add('Vshort_ %s %s step 0' % (Np, Nm))
        else:
            new.add('Vshort_ %s %s 0' % (Np, Nm))

        # Negate current since Vshort is a considered a source.
        Isc = -new.get_I('Vshort_', **kwargs)

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
            return self.norton(Np, Nm)
        except:
            return self.thevenin(Np, Nm)

    def thevenin(self, Np, Nm=None):
        """Return s-domain Thevenin oneport model between nodes Np and Nm.

        If Np is a component name, create model using the component nodes."""

        from .oneport import V, Z

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        Voc = self.Voc(Np, Nm)
        Zoc = self.impedance(Np, Nm)

        return (V(Voc) + Z(Zoc)).simplify()

    def norton(self, Np, Nm=None):
        """Return s-domain Norton model between nodes Np and Nm.

        If Np is a component name, create model using the component nodes."""

        from .oneport import I, Y

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        Isc = self.Isc(Np, Nm)
        Ysc = self.admittance(Np, Nm)

        return (I(Isc) | Y(Ysc)).simplify()

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

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.kill()
        new._add_ground(Nm)
        test = new._add_test_voltage_source(Np, Nm)
        If = new[test].I

        return admittance(If.laplace().sympy)

    def impedance(self, Np, Nm=None):
        """Return driving-point impedance between nodes
        Np and Nm with independent sources killed and initial
        conditions ignored.  Since the result is causal, the frequency
        domain impedance can be found by substituting j * omega for
        s."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.apply_test_current_source(Np, Nm)
        Vf = new.Voc(Np, Nm)
        return impedance(Vf.laplace().sympy)

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
        domain.  See also conductance, reactance, resistance."""
        return self.impedance(Np, Nm).B

    def transfer(self, N1p, N1m, N2p=None, N2m=None):
        """Create s-domain voltage transfer function V2(s) / V1(s) where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed and initial conditions
        are ignored.  Since the result is causal, the frequency response
        can be found by substituting j * omega for s.

        Alternative forms are:
            transfer(N1p, N1m, N2p, N2m)
            transfer(cpt1, cpt2)
            transfer((N1p, N1m), cpt2)
            transfer(cpt1, (N2p, N2m))
        """

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'transfer')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_voltage_source(N1p, N1m)
        V2 = new.Voc(N2p, N2m)
        H = transfer(V2.laplace())
        H.causal = True
        return H

    def voltage_gain(self, N1p, N1m, N2p=None, N2m=None):
        """Create s-domain voltage transfer function V2(s) / V1(s) where:
        V1 is the test voltage applied between N1p and N1m
        V2 is the measured open-circuit voltage between N2p and N2m

        Note, independent sources are killed and initial conditions
        are ignored.  Since the result is causal, the frequency response
        can be found by substituting j * omega for s.

        Alternative forms are:
            voltage_gain(N1p, N1m, N2p, N2m)
            voltage_gain(cpt1, cpt2)
            voltage_gain((N1p, N1m), cpt2)
            voltage_gain(cpt1, (N2p, N2m))
        """

        # This is the same as transfer.
        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'voltage_gain')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_voltage_source(N1p, N1m)
        V2 = new.Voc(N2p, N2m)
        H = transfer(V2.laplace())
        H.causal = True
        return H

    def current_gain(self, N1p, N1m, N2p=None, N2m=None):
        """Create s-domain current transfer function I2(s) / I1(s) where:
        I1 is the test current applied between N1p and N1m
        I2 is the measured short-circuit current flowing from N2m to N2p

        Note, the currents are considered to be flowing into the
        positive nodes as is the convention with two-ports.  Thus the
        input and output currents have opposite directions and so a
        piece of wire has a current gain of -1.

        Note, independent sources are killed and initial conditions
        are ignored.  Since the result is causal, the frequency response
        can be found by substituting j * omega for s.

        Alternative forms are:
            current_gain(N1p, N1m, N2p, N2m)
            current_gain(cpt1, cpt2)
            current_gain((N1p, N1m), cpt2)
            current_gain(cpt1, (N2p, N2m))

        """

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'current_gain')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_current_source(N1p, N1m)
        H = transfer(-new.Isc(N2p, N2m).laplace())
        H.causal = True
        return H

    def transadmittance(self, N1p, N1m, N2p=None, N2m=None):
        """Create s-domain transadmittance (transfer admittance) function
        I2(s) / V1(s) where:
          V1 is the test voltage applied between N1p and N1m
          I2 is the measured short-circuit current flowing from N2m to N2p.

        Note, I2 is considered to be flowing into the positive node as
        is the convention with two-ports.  Thus the transadmittance of
        a series resistor with resistance R is -1 / R.

        Note, independent sources are killed and initial conditions
        are ignored.  Since the result is causal, the frequency response
        can be found by substituting j * omega for s.

        Alternative forms are:
            transadmittance(N1p, N1m, N2p, N2m)
            transadmittance(cpt1, cpt2)
            transadmittance((N1p, N1m), cpt2)
            transadmittance(cpt1, (N2p, N2m))

        """

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'transadmittance')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_voltage_source(N1p, N1m)
        H = admittance(-new.Isc(N2p, N2m).laplace())
        H.causal = True
        return H

    def transimpedance(self, N1p, N1m, N2p=None, N2m=None):
        """Create s-domain transimpedance (transfer impedance) function
        V2(s) / I1(s) where:
          I1 is the test current applied between N1p and N1m
          V2 is the measured open-circuit voltage between N2p and N2m.

        Note, I1 is considered to be flowing into the positive node as
        is the convention with two-ports.

        Note, independent sources are killed and initial conditions
        are ignored.  Since the result is causal, the frequency response
        can be found by substituting j * omega for s.

        Alternative forms are:
            transimpedance(N1p, N1m, N2p, N2m)
            transimpedance(cpt1, cpt2)
            transimpedance((N1p, N1m), cpt2)
            transimpedance(cpt1, (N2p, N2m))
        """

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'transadmittance')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_current_source(N1p, N1m)
        H = impedance(new.Voc(N2p, N2m).laplace())
        H.causal = True
        return H

    def Aparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create A-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Bparams, Gparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """

        from .twoport import AMatrix

        N1p, N1m, N2p, N2m = self._parse_node_args4(
            N1p, N1m, N2p, N2m, 'Aparams')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)
        new = self.kill()
        new._add_ground(N1m)

        try:
            test = new._add_test_voltage_source(N1p, N1m)

            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = new.Voc(N1p, N1m)(s) / new.Voc(N2p, N2m)(s)

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure -I2 with port 2 short-circuit
            A12 = new.Voc(N1p, N1m)(s) / new.Isc(N2p, N2m)(s)

            new.remove(test)

            test = new._add_test_current_source(N1p, N1m)

            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure V2 with port 2 open-circuit
            try:
                A21 = current(0 * s + 1) / new.Voc(N2p, N2m)(s)
            except ValueError:
                # It is likely there is an open-circuit.
                new2 = new.copy()
                new2.add('W %s %s' % (N2p, N2m))
                A21 = -new2[test].I(s) / new2.Voc(N2p, N2m)(s)
                A21 = 0

            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure -I2 with port 2 short-circuit
            A22 = current(0 * s + 1) / new.Isc(N2p, N2m)(s)

            new.remove(test)
            A = AMatrix(((A11, A12), (A21, A22)))
            return A

        except ValueError:
            warn('Cannot create A matrix directly; trying via Z matrix')
            Z = self.Zparams(N1p, N1m, N2p, N2m)
            return Z.Aparams

    def Bparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create B-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Gparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Bparams

    def Gparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create G-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Gparams

    def Hparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create H-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Hparams

    def Sparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create S-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Sparams

    def Tparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create T-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Sparams, Yparams, and Zparams.
        """
        return self.Tparams(N1p, N1m, N2p, N2m).Hparams

    def Yparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create Y-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Sparams, Tparams, and Zparams.
        """
        return self.Zparams(N1p, N1m, N2p, N2m).Yparams

    def Zparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create Z-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Sparams, Tparams, and Yparams.
        """
        from .twoport import ZMatrix

        # TODO, generalise to multiports.

        N1p, N1m, N2p, N2m = self._parse_node_args4(
            N1p, N1m, N2p, N2m, 'Zparams')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)
        new = self.kill()
        new._add_ground(N1m)

        try:
            test = new._add_test_current_source(N1p, N1m)

            # Z11 = V1 / I1 with I2 = 0
            # Apply I1 and measure V1 with port 2 open-circuit
            Z11 = impedance(new.Voc(N1p, N1m)(s))

            # Z21 = V2 / I1 with I2 = 0
            # Apply I1 and measure V2 with port 2 open-circuit
            Z21 = impedance(new.Voc(N2p, N2m)(s))

            new.remove(test)

            test = new._add_test_current_source(N2p, N2m)

            # Z12 = V1 / I2 with I1 = 0
            # Apply I2 and measure V1 with port 1 open-circuit
            Z12 = impedance(new.Voc(N1p, N1m)(s))

            # Z22 = V2 / I2 with I1 = 0
            # Apply I2 and measure V2 with port 1 open-circuit
            Z22 = impedance(new.Voc(N2p, N2m)(s))

            new.remove(test)

            Z = ZMatrix(((Z11, Z12), (Z21, Z22)))
            return Z

        except ValueError as e:
            raise ValueError('Cannot create Z matrix: %s' % e)

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

        new = self.kill()
        new._add_ground(nodes[1])

        try:

            Y = Matrix.zeros(len(ports))

            for col in range(len(ports)):

                for row in range(len(ports)):
                    if row == col:
                        new.add('V%d_ %s %s {DiracDelta(t)}' % (
                            row, ports[row][0], ports[row][1]))
                    else:
                        new.add('V%d_ %s %s 0' %
                                (row, ports[row][0], ports[row][1]))

                for row in range(len(ports)):
                    Y[row, col] = admittance(new.elements['V%d_' % row].I(s))

                for row in range(len(ports)):
                    new.remove('V%d_' % row)
            return Y

        except ValueError as e:
            raise ValueError('Cannot create Y matrix: %s' % e)

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

        new = self.kill()
        new._add_ground(nodes[1])

        try:

            Z = Matrix.zeros(len(ports))

            for col in range(len(ports)):
                new.add('I_ %s %s {DiracDelta(t)}' %
                        (ports[col][0], ports[col][1]))

                for row in range(len(ports)):
                    Z[row, col] = impedance(
                        new.Voc(ports[row][0], ports[row][1])(s))

                new.remove('I_')
            return Z

        except ValueError as e:
            raise ValueError('Cannot create Z matrix: %s' % e)

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
            net = cpt._select(kind)
            new._add(net)
        return new

    def open_circuit(self, cpt):
        """Apply open-circuit in series with component.  Returns name of open
        circuit component."""

        if isinstance(cpt, Cpt):
            cpt = cpt.name
        return self.elements[cpt].open_circuit()

    def short_circuit(self, cpt):
        """Apply short-circuit across component.  Returns name of voltage
        source component used as the short."""

        if isinstance(cpt, Cpt):
            cpt = cpt.name
        return self.elements[cpt].short_circuit()

    def convert_IVP(self, t=0):
        """Remove switches from netlist and convert to an initial value
        problem. `t` is used to determine the state of the switches."""

        times = self.switching_times()
        if times == ():
            warn('Netlist has no switches')
            return self

        if t < times[0]:
            return self.replace_switches(t)

        cct = self
        for m, time in enumerate(times):
            if time > t:
                break
            before = cct.replace_switches_before(time)
            cct = cct.replace_switches(time).initialize(before, time)

        if time != 0:
            warn('Note, the time t is relative to %s' % time)

        return cct

    def _replace_switches(self, t=0, switchnames=None, before=False):
        """Replace specified switches with wire or open circuit
        for time `t`.   If `switchnames` is not specified, all switches
        are replaced."""

        new = self._new()

        for cpt in self._elements.values():
            if switchnames is None or cpt.name in switchnames:
                net = cpt._replace_switch(t, before=before)
            else:
                net = cpt._copy()
            new._add(net)
        return new

    def replace_switches(self, t=0, switchnames=None):
        """Replace specified switches with open-circuit or short-circuit for
        time at or after `t`.  If `switchnames` is not specified, all
        switches are replaced."""

        return self._replace_switches(t, switchnames, before=False)

    def replace_switches_before(self, t=0, switchnames=None):
        """Replace specified switches with open-circuit or short-circuit for
        time just before `t`.  If `switchnames` is not specified, all
        switches are replaced."""

        return self._replace_switches(t, switchnames, before=True)

    def switching_times(self, tmax=1e12):
        """Return sorted list of the times that switches activate prior to
        `tmax`."""

        times = []
        for cpt in self._elements.values():
            if cpt.type.startswith('SW'):
                active_time = float(cpt.args[0])
                if active_time < tmax and active_time not in times:
                    times.append(active_time)
        return sorted(times)

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
                raise ValueError(
                    'Element %s is not a known independent source' % arg)
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
                raise ValueError(
                    'Element %s is not a known independent source' % arg)

        return self._kill(sources)

    def kill_noise(self):
        """Return a new circuit with the independent noise voltage sources and
        noise current sources killed."""

        new = self._new()

        for cpt in self._elements.values():
            if cpt.independent_source and cpt.is_noisy:
                net = cpt._kill()
            else:
                net = cpt._copy()
            new._add(net)
        return new

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
        """Generate circuit graph for this netlist.   This is cached."""

        return self.circuit_graph()

    def circuit_graph(self):
        """Generate circuit graph for this netlist.   This is cached."""

        from .circuitgraph import CircuitGraph

        if hasattr(self, '_cg'):
            return self._cg

        self._cg = CircuitGraph.from_circuit(self)
        return self._cg

    def mesh_analysis(self):
        """Perform mesh analysis for this netlist.   This is cached.

        This is only applicable for circuits with a planar topology.
        """

        from .loopanalysis import LoopAnalysis

        if hasattr(self, '_la'):
            return self._la

        self._la = LoopAnalysis.from_circuit(self)
        return self._la

    def loop_analysis(self):
        """Perform loop analysis for this netlist.   This is cached.

        This is currently an alias for `mesh_analysis()` and so only works
        for circuits with a planar topology.
        """

        return self.mesh_analysis()

    def nodal_analysis(self):
        """Perform nodal analysis for this netlist.   This is cached."""

        from .nodalanalysis import NodalAnalysis

        if hasattr(self, '_na'):
            return self._na

        self._na = NodalAnalysis.from_circuit(self)
        return self._na

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

        parallel_set = self.cg.in_parallel(cpt)
        return parallel_set

    def _in_series_names(self, cpt=None):

        if isinstance(cpt, Cpt):
            cpt = cpt.name

        series_set = self.cg.in_series(cpt)
        return series_set

    def in_parallel(self, cpt=None):
        """Return set of cpts in parallel with specified cpt.  If no cpt
        specified, return list of sets of parallel cpts."""

        if cpt is None:
            return self._in_parallel_all()

        if isinstance(cpt, Cpt):
            cpt = cpt.name

        names = self._in_parallel_names(cpt)
        if len(names) < 2:
            return set()
        names.discard(cpt)
        return list(names)

    def in_series(self, cpt=None):
        """Return set of cpts in series with specified cpt.  If no cpt
        specified, return list of sets of series cpts."""

        if cpt is None:
            return self._in_series_all()

        if isinstance(cpt, Cpt):
            cpt = cpt.name

        names = self._in_series_names(cpt)
        if len(names) < 2:
            return set()
        names.discard(cpt)
        return list(names)

    def check(self):
        """Check if network contains a loop of voltage sources or a cut set of
        current sources."""

        return self.simplify(explain=True, modify=False)

    def twoport(self, N1p, N1m, N2p=None, N2m=None, model='B'):
        """Create s-domain twoport model for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        `model` is `A, `B`, `G`, `H`, `Y`, or `Z`.

        Alternative forms are:
            twoport(N1p, N1m, N2p, N2m)
            twoport(cpt1, cpt2)
            twoport((N1p, N1m), cpt2)
            twoport(cpt1, (N2p, N2m))
        """

        from .twoport import TwoPortAModel, TwoPortBModel, TwoPortGModel
        from .twoport import TwoPortHModel, TwoPortYModel, TwoPortZModel

        N1p, N1m, N2p, N2m = self._parse_node_args4(
            N1p, N1m, N2p, N2m, 'twoport')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        # TODO, generalise for not just s-domain.

        new = self.copy()
        new._add_ground(N1m)

        if model == 'A':
            V1a = new.Voc(N1p, N1m, nowarn=True)(s)
            I1a = new.Isc(N1p, N1m, nowarn=True)(s)
            A = new.Aparams(N1p, N1m, N2p, N2m)
            return TwoPortAModel(A, V1a=V1a, I1a=I1a)
        elif model == 'B':
            V2b = new.Voc(N2p, N2m, nowarn=True)(s)
            I2b = new.Isc(N2p, N2m, nowarn=True)(s)
            A = new.Aparams(N1p, N1m, N2p, N2m)
            return TwoPortBModel(A.Bparams, V2b=V2b, I2b=I2b)
        elif model == 'Z':
            V1 = new.Voc(N1p, N1m, nowarn=True)(s)
            V2 = new.Voc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortZModel(Z, V1z=V1, V2z=V2)
        elif model == 'Y':
            I1 = new.Isc(N1p, N1m, nowarn=True)(s)
            I2 = new.Isc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortYModel(Z.Yparams, I1y=I1, I2y=I2)
        elif model == 'G':
            I1 = new.Isc(N1p, N1m, nowarn=True)(s)
            V2 = new.Voc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortGModel(Z.Gparams, I1g=I1, V2g=V2)
        elif model == 'H':
            V1 = new.Voc(N1p, N1m, nowarn=True)(s)
            I2 = new.Isc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortHModel(Z.Hparams, V1h=V1, I2h=I2)
        else:
            raise ValueError('Model %s unknown, must be B, H, Y, or Z' % model)

    @property
    def sch(self):
        """Generate schematic of subnetlist."""

        from .schematic import Schematic

        if hasattr(self, '_sch'):
            return self._sch

        sch = Schematic(allow_anon=self.allow_anon)

        netlist = self.netlist()
        for net in netlist.split('\n'):
            sch.add(net)

        sch.subnetlists = self.subnetlists
        self._sch = sch
        return sch

    @property
    def sim(self):
        """Generate simulation object."""

        if hasattr(self, '_sim'):
            return self._sim

        self._sim = Simulator(self)
        return self._sim

    def state_space(self, node_voltages=None, branch_currents=None):
        """Generate state-space representation.

        `node_voltages` is a list of node names to use as voltage outputs.
        If `None` use all the unique node names.

        `branch_currents` is a list of component names to use as
        current outputs.  If `None` use all the components.

        Here's an example:
        `cct = Circuit('cct.sch')
        ss = cct.state_space(node_voltages=['1', '3'], branch_currents=['L1', 'L2'])`
        """

        ss = StateSpace.from_circuit(self, node_voltages, branch_currents)
        return ss

    @property
    def ss(self):
        """Generate state-space representation.  See also `state_space()`"""
        return self.state_space()

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

    def pre_initial_model(self):
        """Generate model for determining the pre-initial conditions."""

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

    def _initialize_from_dict(self, values):

        for key, value in values.items():
            if key not in self.reactances:
                warn('Cannot set initial value for %s since not a reactance' % key)
            elif key not in self._elements:
                raise ValueError('Unknown compoment %s' % key)

        new = self._new()

        for cpt in self._elements.values():
            ic = 0
            if cpt.name in values:
                ic = values[cpt.name]
                try:
                    ic = ic.remove_condition()
                except:
                    pass

            net = cpt._initialize(ic)
            new._add(net)
        return new

    def _initialize_from_circuit(self, cct, T=None):

        if T is None:
            raise ValueError('Time T not specified')

        new = self._new()

        for cpt in self._elements.values():
            ic = 0
            if cpt.name in cct.reactances:
                if cpt.type == 'C':
                    ic = cct[cpt.name].v.remove_condition().subs(T)
                else:
                    ic = cct[cpt.name].i.remove_condition().subs(T)

            net = cpt._initialize(ic)
            new._add(net)
        return new

    def initialize(self, cct, T=None):
        """Set the initial values for this netlist based on the values
        computed for netlist `cct` at specified time `T`.

        Alternatively, set the initial values using a dictionary
        of values keyed by the component name.
        """

        if isinstance(cct, dict):
            return self._initialize_from_dict(cct)

        return self._initialize_from_circuit(cct, T)

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
        return self.analysis.causal

    @property
    def is_dc(self):
        """Return True if all independent sources are DC and not an
        initial value problem.  The initial value problem may collapse
        to a DC problem but we cannot prove this yet."""
        return self.analysis.dc

    @property
    def is_ac(self):
        """Return True if all independent sources are AC and not an
        initial value problem."""
        return self.analysis.ac

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
        return self.analysis.has_dc

    @property
    def has_ac(self):
        """Return True if any independent source has an AC component."""
        return self.analysis.has_ac

    @property
    def has_s_transient(self):
        """Return True if any independent source has a transient component defined in s-domain."""
        return self.analysis.has_s

    @property
    def has_transient(self):
        """Return True if any independent source has a transient component."""
        return self.analysis.has_transient

    @property
    def zeroic(self):
        """Return True if the initial conditions for all components are zero."""
        return self.analysis.zeroic

    @property
    def is_IVP(self):
        """Return True for an initial value problem.  This is True if any
        component (L, C) that allows initial conditions has them explicitly
        defined.

        """
        return self.analysis.ivp

    @property
    def is_ivp(self):
        LcapyDeprecationWarning(
            feature="is_ivp",
            useinstead="is_IVP",
            deprecated_since_version="1.7"
        ).warn()
        return self.is_IVP

    @property
    def is_passive(self):
        """Return True for a passive network, i.e., there are no sources."""

        return self.sources == {}

    @property
    def is_switching(self):
        """Return True for a switching circuit."""
        return self.analysis.switching

    @property
    def is_time_domain(self):
        """Return True if can analyse in time domain."""
        return self.analysis.time_domain

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

        return self.analysis.reactances

    @property
    def transformers(self):
        """Return dictionary of transformers."""

        return self.components.transformers

    @property
    def capacitors(self):
        """Return dictionary of capacitors."""

        return self.components.capacitors

    @property
    def inductors(self):
        """Return dictionary of inductors."""

        return self.components.inductors

    @property
    def voltage_sources(self):
        """Return dictionary of voltage_sources."""

        return self.components.voltage_sources

    @property
    def current_sources(self):
        """Return dictionary of current_sources."""

        return self.components.current_sources

    @property
    def ics(self):
        """Return dictionary of components with initial conditions."""

        return self.analysis.ics

    @property
    def independent_sources(self):
        """Return dictionary of independent sources (this does not include
        implicit sources due to initial conditions)."""

        return self.analysis.independent_sources

    @property
    def dependent_sources(self):
        """Return dictionary of dependent sources."""

        return self.analysis.dependent_sources

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
        return self.analysis.control_sources

    @property
    def analysis(self):

        if not hasattr(self, '_analysis'):
            self._analysis = self.analyse()

        return self._analysis

    def analyse(self):

        return Analysis(self, self.components)

    @property
    def components(self):

        if not hasattr(self, '_components'):
            self._components = Components(self)

        return self._components

    def description(self):
        """Return a message describing how circuit is solved."""

        def describe_sources(sources):
            sources_string = ', '.join(sources)
            if len(sources) == 1:
                return 'source %s' % sources_string
            return 'sources %s' % sources_string

        def describe_analysis(method, sources):
            return '%s analysis is used for %s.\n' % (method,
                                                      describe_sources(sources))

        if self.is_switching:
            return '''This has switches and thus is time variant.  Use the convert_IVP(t) method to convert to an initial value problem, specifying the time when to evaluate the switches.'''

        groups = self.independent_source_groups(
            transform=not self.is_time_domain)

        if groups == {}:
            return 'There are no non-zero independent sources so everything is zero.\n'

        if self.is_IVP:
            return 'This has initial conditions for %s so is an initial value '
            'problem solved in the s-domain using Laplace transforms.\n' \
                % ', '.join(self.ics)
            return

        s = ''
        if len(groups) > 1:
            s = 'This is solved using superposition.\n'

        for kind, sources in groups.items():
            if not isinstance(kind, str):
                s += describe_analysis('Phasor', sources)
            elif kind[0] == 'n':
                s += describe_analysis('Noise', sources)
            elif kind == 'dc':
                s += describe_analysis('DC', sources)
            elif kind == 's':
                s += describe_analysis('Laplace', sources)
            elif kind in ('t', 'time'):
                s += describe_analysis('Time-domain', sources)
        return s

    def describe(self):
        """Print a message describing how circuit is solved."""
        print(self.description())

    def Vname(self, name):
        return Vname(name, self.kind)

    def Iname(self, name):
        return Iname(name, self.kind)

    def annotate(self, cpts, *args, **kwargs):
        """Annotate a particular component (or list of components)
        with specified schematic attributes and return new netlist.

        For example:
        `cct.annotate('R1', color='blue')`
        `cct.annotate('R1', 'color=blue, dashed')`
        `cct.annotate(('U1', 'U2'), fill='blue')`

        See also `highlight`."""

        if not isinstance(cpts, (tuple, list, set)):
            cpts = [cpts]

        names = []
        for cpt in cpts:
            if isinstance(cpt, Cpt):
                name = cpt.name
            else:
                name = cpt
            if not self.has(name):
                raise ValueError('Unknown component %s' % name)
            names.append(name)

        new = self._new()

        for cpt in self._elements.values():
            if cpt.name in names:
                new.add(cpt.annotate(*args, **kwargs))
            else:
                new.add(cpt._copy())
        return new

    def highlight(self, cpts, color='blue'):
        """Highlight a particular component (or list of components)
        with specified color and return new netlist.

        See also `annotate`."""

        return self.annotate(cpts, color=color)

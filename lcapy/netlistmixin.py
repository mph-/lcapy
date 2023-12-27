"""This module provides the NetlistMixin class.  This is used for
Netlist and SubNetlist.

Copyright 2020--2023 Michael Hayes, UCECE

"""

from .analysis import Analysis
from .attrdict import AttrDict
from .components import Components
from .current import Iname
from .deprecation import LcapyDeprecationWarning
from .equipotentialnodes import EquipotentialNodes
from .expr import expr
from .mnacpts import Cpt
from .netfile import NetfileMixin
from .netlistnamespace import NetlistNamespace
from .node import Node
from .nodes import Nodes
from .state import state
from .symbols import j, s, omega
from .voltage import Vname
from . import mnacpts
from .cache import lru_cache, cached_property
from collections import OrderedDict
from warnings import warn


class NetlistMixin(object):

    def __init__(self, filename=None, context=None, allow_anon=False):

        self._elements = OrderedDict()
        self.namespaces = {}
        self.nodes = Nodes()
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

    @property
    def analysis(self):

        return self.analyse()

    @cached_property
    def branch_list(self):
        """Determine list of names of branch elements, e.g.,
        ['C1', 'V1', 'R1', 'R2']."""

        branch_list = []
        for key, elt in self.elements.items():
            if elt.type not in ('W', 'O', 'P', 'K', 'XX'):
                branch_list.append(elt.name)
        return branch_list

    @property
    def capacitors(self):
        """Return list of capacitors."""

        return self.components.capacitors

    @property
    def cg(self):
        """Generate circuit graph for this netlist.  This is cached."""

        return self.circuit_graph()

    @property
    def components(self):

        if not hasattr(self, '_components'):
            self._components = Components(self)

        return self._components

    @property
    def control_sources(self):
        """Return list of voltage sources required to specify control
        current for CCVS and CCCS components."""
        return self.analysis.control_sources

    @property
    def current_sources(self):
        """Return list of current_sources."""

        return self.components.current_sources

    @property
    def dependent_sources(self):
        """Return list of dependent sources."""

        return self.analysis.dependent_sources

    @property
    def elements(self):

        if hasattr(self, '_add_elements'):
            if self._elements == {}:
                self._add_elements()

        return self._elements

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
                enodes.add_wire(*elt.node_names)
            elif elt.type.startswith('TL'):
                enodes.add_wire(elt.node_names[1], elt.node_names[3])
            elif elt.type.startswith('TP'):
                enodes.add_wire(elt.node_names[1], elt.node_names[3])
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
        """Return True if any independent source has a transient component defined in Laplace-domain."""
        return self.analysis.has_s

    @property
    def has_transient(self):
        """Return True if any independent source has a transient component."""
        return self.analysis.has_transient

    @property
    def is_causal(self):
        """Return True if all independent sources are causal and not an
        initial value problem (unless all the initial values are zero)."""
        return self.analysis.causal

    @property
    def is_connected(self):
        """Return True if all components are connected."""

        return self.cg.is_connected

    @property
    def is_dc(self):
        """Return True if all independent sources are DC and not an
        initial value problem.  The initial value problem may collapse
        to a DC problem but we cannot prove this yet."""
        return self.analysis.dc

    @property
    def ics(self):
        """Return list of components with initial conditions."""

        return self.analysis.ics

    @property
    def independent_sources(self):
        """Return list of independent sources (this does not include
        implicit sources due to initial conditions)."""

        return self.analysis.independent_sources

    def independent_source_groups(self, transform=False):
        """Return list of source groups.  Each group is a list of
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
            if not elt.is_independent_source:
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
    def inductors(self):
        """Return list of inductors."""

        return self.components.inductors

    @property
    def is_ac(self):
        """Return True if all independent sources are AC and not an
        initial value problem."""
        return self.analysis.ac

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
    def is_superposition(self):
        """Return True if netlist needs to be solved using multiple approaches,
        e.g., ac and dc"""

        return ((self.has_ac and self.has_dc) or
                (self.has_ac and self.has_transient) or
                (self.has_dc and self.has_transient))

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
    def mutual_couplings(self):
        """Return list of mutual couplingss."""

        return self.analysis.mutual_couplings

    @cached_property
    def node_list(self):
        """Determine list of sorted unique node names, e.g.,
        ['0', '1', '2']."""

        # Extract unique nodes.
        node_list = list(self.equipotential_nodes.keys())
        node_list = sorted(node_list)
        # Ensure node '0' is first in the list.
        if '0' in node_list:
            node_list.insert(0, node_list.pop(node_list.index('0')))

        return node_list

    @cached_property
    def node_map(self):
        """Create dictionary mapping node names to the unique
        equipotential node names."""

        enodes = self.equipotential_nodes

        # Create inverted dictionary that maps the node names
        # to the equipotential node names.
        node_map = {}
        for key, nodes in enodes.items():
            for node in nodes:
                node_map[node] = key

        return node_map

    @property
    def reactances(self):
        """Return list of reactances."""

        return self.analysis.reactances

    @property
    def sch(self):
        """Generate schematic of subnetlist."""

        # Don't cache schematic; need to remove splitting of implicit nodes
        # from the draw method otherwise a.draw(); a.draw() fails.

        from .schematic import Schematic

        sch = Schematic(allow_anon=self.allow_anon)

        netlist = self.netlist()
        for net in netlist.split('\n'):
            sch.add(net)

        sch.subnetlists = self.subnetlists

        return sch

    @property
    def sources(self):
        """Return list of all sources (this does not include
        implicit sources due to initial conditions)."""

        return self.dependent_sources + self.independent_sources

    @property
    def transformers(self):
        """Return list of transformers."""

        return self.components.transformers

    @property
    def voltage_sources(self):
        """Return list of voltage_sources."""

        return self.components.voltage_sources

    @property
    def zeroic(self):
        """Return True if the initial conditions for all components are zero."""
        return self.analysis.zeroic

    def _check_nodes(self, *nodes):

        str_nodes = []
        for node in nodes:
            if isinstance(node, int):
                node = '%s' % node
            if node not in self.nodes:
                raise ValueError('Unknown node %s' % node)
            str_nodes.append(node)
        return str_nodes

    def _dummy_node_name(self):
        """Create a dummy node name."""

        return '_' + self._make_anon_node_name()

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

    def _invalidate(self):

        caches = ('nodal_analysis', 'mesh_analysis', 'modified_nodal_analysis',
                  'circuit_graph', '_subs_make', 'analyse')
        for cache in caches:
            getattr(self, cache).cache_clear()

        # Cached properties
        try:
            del self.branch_list
        except:
            pass
        try:
            del self.node_list
        except:
            pass
        try:
            del self.Vdict
        except:
            pass
        try:
            del self.Idict
        except:
            pass
        try:
            del self.node_map
        except:
            pass

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

    def _noisy(self, resistornames, T='T'):

        new = self._new()

        for cpt in self._elements.values():
            if cpt.name in resistornames:
                net = cpt._noisy(T=T)
            else:
                net = cpt._copy()
            new._add(net)
        return new

    def _parse_node_args2(self, Np, Nm=None):

        if Nm is None:
            cpt = self[Np]
            if isinstance(cpt, Node):
                Np, Nm = cpt.name, 0
            else:
                Np, Nm = cpt.node_names[0:2]
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
                if name != 'as_ladder' and arg1.is_voltage_source:
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

    def _potential_combine_names(self):

        names = []
        for name, elt in self.elements.items():
            if elt.type in ('V', 'I', 'R', 'NR', 'C', 'L', 'Y', 'Z'):
                names.append(name)
        return names

    @lru_cache(1)
    def analyse(self):

        return Analysis(self, self.components)

    def across_nodes(self, node1, node2):
        """Return set of components that are connected directly across nodes
        `node1` and `node2`.

        """

        return self.cg.across_nodes(node1, node2)

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

    def check(self):
        """Check if network contains a loop of voltage sources or a cut set of
        current sources."""

        return self.simplify(explain=True, modify=False)

    @lru_cache(1)
    def circuit_graph(self):
        """Generate circuit graph for this netlist.   This is cached."""

        from .circuitgraph import CircuitGraph

        return CircuitGraph.from_circuit(self)

    def copy(self):
        """Create a copy of the netlist"""

        new = self._new()

        for cpt in self._elements.values():
            new._add(cpt._copy())
        return new

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

    def has(self, cpt):
        """Return True if cpt in elements."""

        return cpt in self.elements

    def in_parallel(self, cpt=None):
        """Return set of component names in parallel with specified component
        name.  If no component name specified, return list of sets of
        parallel cpts.

        """

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
        """Return set of component_names in series with specified component
        name.  If no component name specified, return list of sets of
        series cpts.

        """

        if cpt is None:
            return self._in_series_all()

        if isinstance(cpt, Cpt):
            cpt = cpt.name

        names = self._in_series_names(cpt)
        if len(names) < 2:
            return set()
        names.discard(cpt)
        return list(names)

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

    def match(self, pattern):
        """Return list of components names matching regular
        expression pattern."""

        from re import match

        elts = []
        for cpt in self._elements.values():
            if match(pattern, cpt.name):
                elts.append(cpt.name)
        return elts

    def netlist(self):
        """Return the current netlist."""

        return '\n'.join([str(cpt) for cpt in self._elements.values()])

    def pdb(self):
        """Enter the python debugger."""

        import pdb
        pdb.set_trace()
        return self

    def subs(self, *args):
        """Substitute symbols in the netlist.

        `args` is either:
        - two arguments, e.g. foo.subs(old, new)
        - one dict or set argument whose key/value items correspond to
          old/new pairs.

        For example, `b = a.subs({'R1': 1e3, 'R2': 9e3})`

        Note, this does not substitute component names.
        """

        if len(args) > 2:
            raise ValueError('Too many args')
        if len(args) == 2:
            subs_dict = {args[0]: args[1]}
        else:
            subs_dict = args[0]

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
        self.kind = 'dc'

        for cpt in self._elements.values():
            net = cpt._r_model()
            new._add(net)
        return new

    def s_model(self, kind='s'):
        """"Create Laplace-domain model."""

        new = self._new()
        new.kind = kind

        for cpt in self._elements.values():
            net = cpt._s_model(kind)
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

    def unconnected_nodes(self):
        """Return list of node names that are not connected."""

        return [node.name for node in self.nodes.values() if node.count <= 1]

    def Vname(self, name):
        return Vname(name, self.kind)

    def Iname(self, name):
        return Iname(name, self.kind)

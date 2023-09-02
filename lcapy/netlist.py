"""This module provides the Netlist class.  It could be rolled into
the Circuit class.

Copyright 2014--2023 Michael Hayes, UCECE

"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from .admittance import admittance
from .impedance import impedance
from .config import solver_method
from .current import Iname, current
from .deprecation import LcapyDeprecationWarning
from .expr import Expr, expr
from .matrix import Matrix
from .mna import Nodedict, Branchdict
from .mnacpts import Cpt
from .netlistmixin import NetlistMixin
from .netlistsimplifymixin import NetlistSimplifyMixin
from .simulator import Simulator
from .statespace import StateSpace
from .subnetlist import SubNetlist
from .superpositionvoltage import SuperpositionVoltage
from .superpositioncurrent import SuperpositionCurrent
from .symbols import j, s, t, omega
from .transfer import transfer
from .transformdomains import TransformDomains
from .utils import isiterable
from copy import copy
from warnings import warn


class Netlist(NetlistMixin, NetlistSimplifyMixin):
    """This class handles a generic netlist with multiple sources.
    During analysis, subnetlists are created for each source kind (dc,
    ac, transient, etc).  Since linearity is assumed, superposition is
    employed.

    """

    def __init__(self, filename=None, context=None, allow_anon=False):

        super(Netlist, self).__init__(filename, context, allow_anon=allow_anon)
        self._invalidate()
        self.kind = 'super'
        self.solver_method = solver_method

    def _groups(self):
        """Return dictionary of source groups keyed by domain.

        If the netlist is for an initial value problem, all the
        sources are in a single group called 'ivp'.  Any noise sources
        are ignored.

        If the netlist can be solved in the time-domain (i.e., there
        are no reactive components), all the sources are in single
        group called 'time'.

        Otherwise, the sources are decomposed and then grouped into
        the 'dc', 's', 'n*', and omega categories.  Note, a source can
        appear in multiple groups, for example, a source with voltage
        3 + u(t) will appear in both the 'dc' and 's' groups.

        """

        if self.is_IVP:

            def namelist(elements):
                return ', '.join([elt for elt in elements])

            if self.missing_ic != {}:
                warn('Missing initial conditions for %s' %
                     namelist(self.missing_ic))

            groups = self.independent_source_groups()
            newgroups = {'ivp': []}
            for key, sources in groups.items():
                if isinstance(key, str) and key[0] == 'n':
                    warn('Ignoring noise source %s for initial value problem' % sources)
                else:
                    newgroups['ivp'] += sources
            return newgroups

        elif self.is_time_domain:

            groups = self.independent_source_groups()
            newgroups = {'time': []}
            for key, sources in groups.items():
                if isinstance(key, str) and key[0] == 'n':
                    newgroups[key] = sources
                else:
                    newgroups['time'] += sources
            return newgroups

        else:
            return self.independent_source_groups(transform=True)

    def _subs_make(self, nowarn=False):

        if hasattr(self, '_sub'):
            return self._sub

        cct = self.expand()
        groups = cct._groups()
        sub = TransformDomains()

        for kind, sources in groups.items():
            sub[kind] = SubNetlist(cct, kind)

        if sub == {} and not nowarn:
            warn('Netlist has no sources')

        self._sub = sub
        return sub

    @property
    def cpts(self):
        """Return list of component names."""

        # Perhaps should prune wires, open-circuits, etc. ?
        return list(self._elements.keys())

    @property
    def is_connected(self):
        """Return True if all components are connected."""

        return self.cg.is_connected

    @property
    def kinds(self):
        """Return list of transform domain kinds required to analyse the netlist."""
        return list(self.sub.keys())

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
    def sim(self):
        """Generate simulation object."""

        if hasattr(self, '_sim'):
            return self._sim

        self._sim = Simulator(self)
        return self._sim

    @property
    def ss(self):
        """Generate state-space representation.  See also `state_space()`"""
        return self.state_space()

    @property
    def sub(self):
        """Return dictionary of subnetlists keyed by transform domain kind.
        Note, the subnetlists are not created until a specific one is
        selected.

        """
        return self._subs_make()

    @property
    def subcircuits(self):
        """Return dictionary of subnetlists keyed by transform domain kind.
        Note, the subnetlists are not created until a specific one is
        selected.  The subcircuit keys are :
        'ivp' for an initial value problem solved using Laplace methods,
        's' for transient analysis using Laplace methods,
        'dc' for DC analysis,
        'time' for time-domain analysis when there are no reactive
        components,
        'n*' for noise-analysis (there is a subcircuit for
        each independent noise source), and
        omega (where omega is a number of expression specifying the angular
        frequency) for phasor analysis.

        """
        return self.sub

    @property
    def super_nodes(self):
        """Super nodes are nodes linked by voltage sources."""

        snodes = []

        for elt in self.elements.values():
            if not elt.is_voltage_source:
                continue
            snodes.append(elt.node_names[0:2])

        from .utils import merge_common

        # Merge supernodes to create super-supernodes...
        return list(merge_common(snodes))

    @property
    def symbols(self):
        """Return dictionary of symbols defined in the netlist."""

        return self.context.symbols

    @property
    def undefined_symbols(self):

        symbols = []
        for k, elt in self.elements.items():
            cpt = elt.cpt
            for arg in cpt.args:
                symbols += expr(arg).symbols

        return symbols

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
                    result[node] = SuperpositionCurrent()
                result[node].add(value)
        self._Idict = result
        return result

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
                    result[node] = SuperpositionVoltage()
                result[node].add(value)

        self._Vdict = result
        return result

    def _add_ground(self, node=None):

        if '0' in self.nodes:
            return None

        if node is None:
            node = list(self.node_map)[0]
            warn('Ground node not specified: using node ' + node)

        return self.add('W %s 0' % node)

    def _add_test_voltage_source(self, Np, Nm):

        self._add('V? %s %s {DiracDelta(t)}' % (Np, Nm))
        return self.last_added()

    def _add_test_current_source(self, Np, Nm):

        self._add('I? %s %s {DiracDelta(t)}' % (Np, Nm))
        return self.last_added()

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

        self._namespace_add(cpt.namespace)

    def get_I(self, name, nowarn=False):
        """Current through component (time-domain)"""

        self._add_ground()
        subs = self._subs_make(nowarn=nowarn)

        result = SuperpositionCurrent()
        for sub in subs.values():
            I = sub.get_I(name)
            result.add(I)
        result = result
        return result

    def get_i(self, name):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def _get_Vd(self, Np, Nm=None, nowarn=False):
        """This does not check nodes."""

        self._add_ground()
        subs = self._subs_make(nowarn=nowarn)

        result = SuperpositionVoltage()
        for sub in subs.values():
            Vd = sub.get_Vd(Np, Nm)
            result.add(Vd)
        result = result.canonical()
        return result

    def get_Vd(self, Np, Nm=None, **kwargs):
        """Voltage drop between nodes (time-domain)"""

        Np, Nm = self._check_nodes(Np, Nm)
        return self._get_Vd(Np, Nm, **kwargs)

    def get_vd(self, Np, Nm=None):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(Np, Nm).time()

    def ac(self):
        """Return netlist for ac components of independent sources
        for angular frequency omega.

        See also: dc, transient, laplace.
        """
        # Could look at all the ac frequencies and if there is only
        # one use that?  If have multiple ac frequencies should issue
        # warning.
        return self.select(omega)

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

        cct = self.copy()
        groundwire = cct._add_ground()
        if cpts is None:
            cpts = []
            for elt in cct._elements.values():
                if (elt.is_resistor or elt.is_capacitor or
                        elt.is_inductor or elt.is_voltage_source):
                    cpts.append(elt.name)

        label = ('f' if flow else 'i') + pos

        if domainvar is None:
            domainvar = t

        new = self._new()
        for cpt in cct._elements.values():
            net = cpt._copy()
            if cpt.name in cpts:
                I = cpt.I(domainvar)
                if evalf:
                    I = I.evalf(num_digits)
                net += ', ' if ';' in net else '; '
                net += ', %s={$%s$}' % (label, I.latex_with_units(
                    eng_format=eng_format, evalf=evalf, num_digits=num_digits, show_units=show_units))
            new.add(net)
        if groundwire is not None:
            new.remove(groundwire.name)
        return new

    def annotate_voltages(self, cpts=None, domainvar=None,
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

        cct = self.copy()
        groundwire = cct._add_ground()
        if cpts is None:
            cpts = []
            for elt in cct._elements.values():
                if (elt.is_resistor or elt.is_capacitor or
                        elt.is_inductor or elt.is_current_source):
                    cpts.append(elt.name)

        if domainvar is None:
            domainvar = t

        new = self._new()
        for cpt in cct._elements.values():
            net = cpt._copy()
            if cpt.name in cpts:
                V = cpt.V(domainvar)
                if evalf:
                    V = V.evalf(num_digits)
                net += ', ' if ';' in net else '; '
                net += 'v%s={$%s$}' % (pos, V.latex_with_units(eng_format=eng_format,
                                       evalf=evalf, num_digits=num_digits, show_units=show_units))
            new.add(net)
        if groundwire is not None:
            new.remove(groundwire.name)
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

        cct = self.copy()
        groundwire = cct._add_ground()
        if nodes is None:
            nodes = cct.node_list
        elif not isiterable(nodes):
            nodes = (nodes, )

        if domainvar is None:
            domainvar = t

        new = self.copy()
        for node in nodes:
            v = cct[node].V(domainvar)
            if evalf:
                v = v.evalf(num_digits)

            vstr = '%s' % v.latex_with_units(eng_format=eng_format, evalf=evalf,
                                             num_digits=num_digits,
                                             show_units=show_units)

            if label_voltages:
                vstr = 'V_{%s}=' % node + vstr

            new.add('A%s %s; l={%s}, anchor=%s' % (node, node, vstr, anchor))
        if groundwire is not None:
            new.remove(groundwire.name)
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

    def dc(self):
        """Return netlist for dc components of independent sources.

        See also: ac, transient, laplace.
        """
        return self.select('dc')

    def describe(self):
        """Print a message describing how circuit is solved."""
        print(self.description())

    def description(self):
        """Return a message describing how circuit is solved."""

        def describe_sources(sources, omega=None):
            sources_string = ', '.join(sources)
            if len(sources) == 1:
                s = 'source %s' % sources_string
            else:
                s = 'sources %s' % sources_string

            if omega is not None:
                s += ' at angular frequency ' + str(omega)
            return s

        def describe_analysis(method, sources, omega=None):
            return '%s analysis is used for %s.\n' % (method,
                                                      describe_sources(sources, omega))

        if self.is_switching:
            return 'This has switches and thus is time variant.  Use the convert_IVP(t) method to convert to an initial value problem, specifying the time when to evaluate the switches.'

        groups = self.independent_source_groups(
            transform=not self.is_time_domain)

        if groups == {}:
            return 'There are no non-zero independent sources so everything is zero.\n'

        if self.is_IVP:
            return 'This has initial conditions for %s so is an initial value problem solved in the Laplace-domain using Laplace transforms.\n' \
                % ', '.join(self.ics)
            return

        s = ''
        if len(groups) > 1:
            s = 'This is solved using superposition.\n'

        for kind, sources in groups.items():
            if not isinstance(kind, str):
                s += describe_analysis('Phasor', sources, kind)
            elif kind[0] == 'n':
                s += describe_analysis('Noise', sources)
            elif kind == 'dc':
                s += describe_analysis('DC', sources)
            elif kind == 's':
                s += describe_analysis('Laplace', sources)
            elif kind in ('t', 'time'):
                s += describe_analysis('Time-domain', sources)
        return s

    def expand(self):
        """Expand the netlist, replacing complicated components with simpler
        components."""

        new = self._new()

        for cpt in self._elements.values():
            new._add(cpt._expand())
        return new

    def laplace(self):
        """Return netlist for Laplace representations of independent
        source values.

        See also: dc, ac, transient.

        """
        return self.select('laplace')

    def transient(self):
        """Return netlist for transient components of independent
        sources.  Note, unlike the similar laplace method, dc and ac
        components are ignored.

        See also: dc, ac, laplace.

        """
        return self.select('s')

    def _new(self):

        # TODO.  Copy or share?
        context = self.context
        return self.__class__(context=context)

    def remove(self, name):
        """Remove specified element by name or elements specified in list."""

        if name is None:
            return self

        if isinstance(name, (list, tuple)):
            for name1 in name:
                self.remove(name1)
            return self

        if name not in self._elements:
            raise ValueError('Unknown component: ' + name)

        self._invalidate()

        cpt = self._elements[name]
        for node in cpt.nodes:
            node.remove(cpt)

        self._elements.pop(name, None)
        return self

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

    def oneport(self, Np, Nm=None):
        """Return oneport object between nodes Np and Nm.  This might be a
        Thevenin network, a Norton network, or a single component.

        If Np is a component name, create model using the component nodes."""

        try:
            return self.norton(Np, Nm)
        except:
            return self.thevenin(Np, Nm)

    def thevenin(self, Np, Nm=None):
        """Return Laplace-domain Thevenin oneport model between nodes Np and Nm.

        If Np is a component name, create model using the component nodes."""

        from .oneport import V, Z

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        Voc = self.Voc(Np, Nm)
        Zoc = self.impedance(Np, Nm)

        return (V(Voc) + Z(Zoc)).simplify()

    def nodal_analysis(self):
        """Perform nodal analysis for this netlist.   This is cached."""

        from .nodalanalysis import NodalAnalysis

        if hasattr(self, '_na'):
            return self._na

        self._na = NodalAnalysis.from_circuit(self)
        return self._na

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

    def noise_model(self, T='T'):
        """"Create noise model where resistors are converted into a series
        combination of an ideal resistor and a noise voltage
        source."""

        return self.noisy(T=T)

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

    def norton(self, Np, Nm=None):
        """Return Laplace-domain Norton model between nodes Np and Nm.

        If Np is a component name, create model using the component nodes."""

        from .oneport import I, Y

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        Isc = self.Isc(Np, Nm)
        Ysc = self.admittance(Np, Nm)

        return (I(Isc) | Y(Ysc)).simplify()

    def admittance(self, Np, Nm=None):
        """Return driving-point Laplace-domain admittance between nodes
        Np and Nm with independent sources killed and initial
        conditions ignored."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.kill()
        new._add_ground(Nm)
        test = new._add_test_voltage_source(Np, Nm)
        If = new[test].I

        return admittance(If.laplace().sympy)

    def impedance(self, Np, Nm=None):
        """Return driving-point Laplace-domain impedance between nodes
        Np and Nm with independent sources killed and initial
        conditions ignored."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.apply_test_current_source(Np, Nm)
        Vf = new.Voc(Np, Nm)
        return impedance(Vf.laplace().sympy)

    def initialize(self, cct, T=None):
        """Set the initial values for this netlist based on the values
        computed for netlist `cct` at specified time `T`.

        Alternatively, set the initial values using a dictionary
        of values keyed by the component name.
        """

        if isinstance(cct, dict):
            return self._initialize_from_dict(cct)

        return self._initialize_from_circuit(cct, T)

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

    def ladder(self, N1p, N1m, N2p=None, N2m=None):
        """Return two-port unbalanced ladder network or `None` if the netlist
        does not have a ladder topology between the specified nodes.

        The input port is defined by the nodes `N1p` and `N1m`.
        The output port is defined by the nodes `N2p` and `N2m`.

        The nodes `N1p` and `N1m` must be the same.

        Alternative forms are:
            ladder(N1p, N1m, N2p, N2m)
            ladder(cpt1, cpt2)
            ladder((N1p, N1m), cpt2)
            ladder(cpt1, (N2p, N2m))

        """

        from lcapy.laddernetworkmaker import LadderNetworkMaker

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'ladder')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        lm = LadderNetworkMaker(self)

        tp = lm.make(N1p, N1m, N2p, N2m)
        return tp

    def loop_analysis(self):
        """Perform loop analysis for this netlist.   This is cached.

        This is currently an alias for `mesh_analysis()` and so only works
        for circuits with a planar topology.
        """

        return self.mesh_analysis()

    def mesh_analysis(self):
        """Perform mesh analysis for this netlist.   This is cached.

        This is only applicable for circuits with a planar topology.
        """

        from .loopanalysis import LoopAnalysis

        if hasattr(self, '_la'):
            return self._la

        self._la = LoopAnalysis.from_circuit(self)
        return self._la

    def reactance(self, Np, Nm=None):
        """Return reactance between nodes Np and Nm with independent
        sources killed.  The result is in the AC (omega) domain.
        See also conductance, resistance, susceptance."""
        return self.impedance(Np, Nm).X

    def resistance(self, Np, Nm=None):
        """Return resistance between nodes Np and Nm with independent
        sources killed.  The result is in the AC (omega) domain.
        See also conductance, reactance, susceptance."""
        return self.impedance(Np, Nm).R

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

    def transfer(self, N1p, N1m, N2p=None, N2m=None):
        """Create Laplace-domain voltage transfer function V2(s) / V1(s) where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed and initial conditions
        are ignored.

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

    def twoport(self, N1p, N1m, N2p=None, N2m=None, model='B'):
        """Create Laplace-domain twoport model for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
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

        # TODO, generalise for not just Laplace-domain.

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

    def voltage_gain(self, N1p, N1m, N2p=None, N2m=None):
        """Create Laplace-domain voltage transfer function V2(s) / V1(s) where:
        V1 is the test voltage applied between N1p and N1m
        V2 is the measured open-circuit voltage between N2p and N2m

        Note, independent sources are killed and initial conditions
        are ignored.

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
        """Create Laplace-domain current transfer function I2(s) / I1(s) where:
        I1 is the test current applied between N1p and N1m
        I2 is the measured short-circuit current flowing from N2m to N2p

        Note, the currents are considered to be flowing into the
        positive nodes as is the convention with two-ports.  Thus the
        input and output currents have opposite directions and so a
        piece of wire has a current gain of -1.

        Note, independent sources are killed and initial conditions
        are ignored.

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
        """Create Laplace-domain transadmittance (transfer admittance) function
        I2(s) / V1(s) where:
          V1 is the test voltage applied between N1p and N1m
          I2 is the measured short-circuit current flowing from N2m to N2p.

        Note, I2 is considered to be flowing into the positive node as
        is the convention with two-ports.  Thus the transadmittance of
        a series resistor with resistance R is -1 / R.

        Note, independent sources are killed and initial conditions
        are ignored.

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
        """Create Laplace-domain transimpedance (transfer impedance) function
        V2(s) / I1(s) where:
          I1 is the test current applied between N1p and N1m
          V2 is the measured open-circuit voltage between N2p and N2m.

        Note, I1 is considered to be flowing into the positive node as
        is the convention with two-ports.

        Note, independent sources are killed and initial conditions
        are ignored.

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

    def Voc(self, Np, Nm=None, **kwargs):
        """Return open-circuit transform-domain voltage between nodes Np and
        Nm."""

        return self._get_Vd(Np, Nm, **kwargs)

    def voc(self, Np, Nm=None):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).time()

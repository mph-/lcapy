"""This module provides the Netlist class.  It could be rolled into
the Circuit class.

Copyright 2014--2023 Michael Hayes, UCECE

"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from .admittance import admittance
from .config import solver_method
from .current import Iname, current
from .deprecation import LcapyDeprecationWarning
from .expr import Expr, expr, ExprList
from .mna import Nodedict, Branchdict
from .mnacpts import Cpt
from .netlistmixin import NetlistMixin
from .netlistopsmixin import NetlistOpsMixin
from .netlistsimplifymixin import NetlistSimplifyMixin
from .simulator import Simulator
from .subnetlist import SubNetlist
from .superpositionvoltage import SuperpositionVoltage
from .superpositioncurrent import SuperpositionCurrent
from .symbols import j, s, t, omega, omega0
from .transformdomains import TransformDomains
from .utils import isiterable
from .cache import lru_cache, cached_property
from copy import copy
from warnings import warn


class Netlist(NetlistOpsMixin, NetlistMixin, NetlistSimplifyMixin):
    """This class handles a generic netlist with multiple sources.
    During analysis, subnetlists are created for each source kind (dc,
    ac, transient, etc).  Since linearity is assumed, superposition is
    employed.

    """

    def __init__(self, filename=None, context=None,
                 allow_anon=False, kind='super'):

        super(Netlist, self).__init__(filename, context, allow_anon=allow_anon)
        self._invalidate()
        self.kind = kind
        self.solver_method = solver_method

    @property
    def _time_kind(self):

        return self.kind in ('super', 'time', 't')

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

    @lru_cache(1)
    def _subs_make(self, nowarn=False):

        cct = self.expand()
        groups = cct._groups()
        sub = TransformDomains()

        for kind, sources in groups.items():
            sub[kind] = SubNetlist(cct, kind)

        if sub == {} and not nowarn:
            warn('Netlist has no sources')

        return sub

    @property
    def cpts(self):
        """Return list of component names."""

        # Perhaps should prune wires, open-circuits, etc. ?
        return list(self._elements.keys())

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

    @cached_property
    def Idict(self):
        """Return dictionary of branch currents for each transform domain"""

        result = Branchdict()
        for sub in self.sub.values():
            for node, value in sub.Idict.items():
                if node not in result:
                    result[node] = SuperpositionCurrent()
                result[node].add(value)
        return result

    @cached_property
    def Vdict(self):
        """Return dictionary of node voltages for each transform domain"""

        result = Nodedict()
        for sub in self.sub.values():
            for node, value in sub.Vdict.items():
                if node not in result:
                    result[node] = SuperpositionVoltage()
                result[node].add(value)

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

    def ac_omega_list(self):
        """Return list of the angular frequencies of the sources in
        the netlist."""

        omega_list = []
        for group in self.independent_source_groups(True).keys():
            if group in ('dc', 's'):
                continue
            if isinstance(group, str) and group[0] == 'n':
                continue
            omega_list.append(group)
        return omega_list

    def ac(self, omega=None):
        """Return netlist for ac components of independent sources
        for angular frequency `omega`.  If `omega` is undefined,
        the angular frequency `omega0` is used.

        See also: dc, transient, laplace.
        """

        if omega is None:

            omega_list = self.ac_omega_list()
            if len(omega_list) > 1:
                omega = omega_list[0]
                warn('Netlist has multiple AC frequencies: %s, using %s' %
                     (omega_list, omega))
            elif len(omega_list) == 1:
                omega = omega_list[0]
            else:
                # Dummy, it could be anything.
                omega = 1

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
        the current is a number, e.g., 100\, mV instead of 0.1\, V

        `evalf` (default True) if True prints floating point
        numbers as decimals otherwise they are shown as rationals

        `show_units` (default True) if True applies the units(e.g.,
        V for volts)

        `pos` specifies where to position the labels(see docs)
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
        the voltage is a number, e.g., 100\, mV instead of 0.1\, V

        `evalf` (default True) if True prints floating point
        numbers as decimals otherwise they are shown as rationals

        `show_units` (default True) if True applies the units(e.g.,
        V for volts)

        `pos` specifies where to position the labels(see docs)
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
        annotation with V1 = for node 1, etc.

        `eng_format` (default True) if True use engineering format if
        the voltage is a number, e.g., 100\, mV instead of 0.1\, V

        `evalf` (default True) if True prints floating point
        numbers as decimals otherwise they are shown as rationals

        `num_digits` (default 3) specfies the number of digits to print
        for floating point numbers

        `show_units` (default True) if True applies the units(e.g.,
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
        """This copies the netlist, removes any voltage sources across the
        specified nodes, kills all the remaining sources, and applies
        a Dirac delta test voltage source across the specified nodes.
        If the netlist is not connected to ground, the negative
        specified node is connected to ground.  The new netlist is
        returned.

        """

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.copy()
        for name in self.cg.across_nodes(Np, Nm):
            if self[name].is_voltage_source:
                warn('Removing voltage source %s across input nodes %s, %s' %
                     (name, Np, Nm))
                new.remove(name)

        new = new.kill()
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
        return self.__class__(context=context, kind=self.kind)

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

    def matrix_equations(self, form='default', invert=False):
        """Return modified nodal analysis equations in matrix form.

        Forms can be:
         'default'
         'A y = b'
         'b = A y'
         'Ainv b = y'
         'y = Ainv b'

        If `invert` is True, evaluate the matrix inverse."""

        cct = self.expand()
        mna = cct.modified_nodal_analysis()
        return mna.matrix_equations(form, invert)

    @lru_cache(1)
    def modified_nodal_analysis(self):
        """Perform modified nodal analysis for this netlist.  This is
        cached."""

        if self.kind in ('time', 'super'):
            raise ValueError(
                'Cannot put time domain equations into matrix form.  '
                'Convert to dc, ac, or laplace domain first.')

        return SubNetlist(self, self.kind).mna

    @lru_cache(1)
    def nodal_analysis(self, node_prefix=''):
        """Perform nodal analysis for this netlist.   This is cached.

        Use `node_prefix` to distingish between nodal voltages and
        voltage source voltages."""

        from .nodalanalysis import NodalAnalysis

        return NodalAnalysis.from_circuit(self, node_prefix=node_prefix)

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

    def initialize(self, cct, T=None):
        """Set the initial values for this netlist based on the values
        computed for netlist `cct` at specified time `T`.

        Alternatively, set the initial values using a dictionary
        of values keyed by the component name.
        """

        if isinstance(cct, dict):
            return self._initialize_from_dict(cct)

        return self._initialize_from_circuit(cct, T)

    def loop_analysis(self):
        """Perform loop analysis for this netlist.   This is cached.

        This is currently an alias for `mesh_analysis()` and so only works
        for circuits with a planar topology.
        """

        return self.mesh_analysis()

    @lru_cache(1)
    def mesh_analysis(self):
        """Perform mesh analysis for this netlist.   This is cached.

        This is only applicable for circuits with a planar topology.
        """

        from .loopanalysis import LoopAnalysis

        return LoopAnalysis.from_circuit(self)

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
        new.kind = kind

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

    def branch_current_names(self):
        """Return `ExprList` of branch current names of the form i_cptname
        for the time-domain and of the form I_cptname othwerwise."""

        branch_list = self.branch_list

        prefix = 'I'
        if self._time_kind:
            prefix = 'i'

        current_name_list = ExprList()
        for branch in branch_list:
            current_name_list.append(prefix + branch)

        return current_name_list

    def branch_currents(self):
        """Return `ExprList` of branch currents.  Each element is a
        SuperpositionVoltage object.

        If you want a vector of time-domain expressions use
        `Vector(cct.branch_currents()(t))`.
        """

        return ExprList([self[b].I for b in self.branch_list])

    def branch_voltage_names(self):
        """Return `ExprList` of branch voltage names of the form v_cptname
        for the time-domain and of the form V_cptname othwerwise."""

        branch_list = self.branch_list

        prefix = 'V'
        if self._time_kind:
            prefix = 'v'

        voltage_name_list = ExprList()
        for branch in branch_list:
            voltage_name_list.append(prefix + branch)

        return voltage_name_list

    def branch_voltages(self):
        """Return `ExprList` of branch voltages.  Each element is a
        SuperpositionVoltage object.

        If you want a vector of time-domain expressions use
        `Vector(cct.branch_voltages()(t))`.
        """

        return ExprList([self[b].V for b in self.branch_list])

    def evidence_matrix(self):
        """Return the evidence matrix.  This has a size NxM where
        N is the number of nodes and M is the number of branches.

        An element is 1 for currents entering a node from a branch, -1
        for currents leaving a node from a branch, and zero otherwise.

        Note, each column sums to zero.

        """

        from .matrix import Matrix
        from sympy import zeros

        cct = self

        branch_list = cct.branch_list
        current_list = self.branch_currents

        nodes = list(cct.equipotential_nodes)
        node_map = cct.node_map

        mat = Matrix(zeros(len(nodes), len(branch_list)))

        for m, cpt_name in enumerate(branch_list):
            cpt = cct[cpt_name]
            node1, node2 = cpt.nodes[0:2]
            n1 = nodes.index(node_map[node1.name])
            n2 = nodes.index(node_map[node2.name])

            mat[n1, m] = 1
            mat[n2, m] = -1

        return mat

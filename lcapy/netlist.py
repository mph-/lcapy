"""This module provides the Netlist class.  It could be rolled into
the Circuit class.

Copyright 2014--2023 Michael Hayes, UCECE

"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

# TODO: This needs an overhaul to untangle the classes.

from __future__ import division
from .superpositionvoltage import SuperpositionVoltage
from .superpositioncurrent import SuperpositionCurrent
from .current import Iname
from .netlistmixin import NetlistMixin
from .netlistsimplifymixin import NetlistSimplifyMixin
from .expr import Expr, expr
from .subnetlist import SubNetlist
from .mna import Nodedict, Branchdict
from .symbols import omega
from copy import copy
from warnings import warn


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
        sub = Transformdomains()

        for kind, sources in groups.items():
            sub[kind] = SubNetlist(cct, kind)

        if sub == {} and not nowarn:
            warn('Netlist has no sources')

        self._sub = sub
        return sub

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
    def kinds(self):
        """Return list of transform domain kinds required to analyse the netlist."""
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
                    result[node] = SuperpositionVoltage()
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
                    result[node] = SuperpositionCurrent()
                result[node].add(value)
        self._Idict = result
        return result

    def _add_ground(self, node):

        if '0' not in self.nodes:
            self.add('W %s 0' % node)

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

    def dc(self):
        """Returnnetlist for dc components of independent sources.

        See also, ac, transient, laplace.
        """
        return self.select('dc')

    def ac(self):
        """Returnnetlist for ac components of independent sources
        for angular frequency omega.

        See also, dc, transient, laplace.
        """
        # Could look at all the ac frequencies and if there is only
        # one use that?  If have multiple ac frequencies should issue
        # warning.
        return self.select(omega)

    def transient(self):
        """Returnnetlist for transient components of independent
        sources.  Note, unlike the similar laplace method, dc and ac
        components are ignored.

        See also, dc, ac, laplace.

        """
        return self.select('s')

    def laplace(self):
        """Returnnetlist for Laplace representations of independent
        source values.

        See also, dc, ac, transient.

        """
        return self.select('laplace')

    @property
    def undefined_symbols(self):

        symbols = []
        for k, elt in self.elements.items():
            cpt = elt.cpt
            for arg in cpt.args:
                symbols += expr(arg).symbols

        return symbols

    def expand(self):
        """Expand the netlist, replacing complicated components with simpler
        components."""

        new = self._new()

        for cpt in self._elements.values():
            new._add(cpt._expand())
        return new

    def _new(self):

        # TODO.  Copy or share?
        context = self.context
        return self.__class__(context=context)

    def remove(self, name):
        """Remove specified element or elements specified in list."""

        if isinstance(name, (list, tuple)):
            for name1 in name:
                self.remove(name1)
            return self

        self._invalidate()

        if name not in self._elements:
            raise ValueError('Unknown component: ' + name)

        cpt = self._elements[name]
        for node in cpt.nodes:
            node.remove(cpt)

        self._elements.pop(name, None)
        return self

"""This module provides support for the common aspects of Circuit and
Network classes.

Copyright 2014--2021 Michael Hayes, UCECE

"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

# TODO: This needs an overhaul to untangle the classes.

from __future__ import division
from .superpositionvoltage import SuperpositionVoltage
from .superpositioncurrent import SuperpositionCurrent
from .current import Iname
from .schematic import Schematic
from .netlistmixin import NetlistMixin
from .netfile import NetfileMixin
from .expr import Expr, expr
from .subnetlist import SubNetlist
from .mna import MNAMixin, Nodedict, Branchdict
from .symbols import omega
from copy import copy


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
                     '_node_map', '_ss', '_node_list', '_branch_list', '_cg'):
            try:
                delattr(self, attr)
            except:
                pass

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
            
        if self.is_ivp:
            
            def namelist(elements):
                return ', '.join([elt for elt in elements])

            if self.missing_ic != {}:
                print('Warning: missing initial conditions for %s' %
                      namelist(self.missing_ic))

            groups = self.independent_source_groups()                       
            newgroups = {'ivp' : []}
            for key, sources in groups.items():
                if isinstance(key, str) and key[0] == 'n':
                    print('Warning: ignoring noise source %'
                          ' for initial value problem' % sources)
                else:
                    newgroups['ivp'] += sources
            return newgroups

        elif self.is_time_domain:
            
            groups = self.independent_source_groups()            
            newgroups = {'time' : []}
            for key, sources in groups.items():
                if isinstance(key, str) and key[0] == 'n':
                    newgroups[key] = sources
                else:
                    newgroups['time'] += sources
            return newgroups

        else:
            return self.independent_source_groups(transform=True)        

        
    def _sub_make(self):

        groups = self._groups()
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

    def get_I(self, name):
        """Current through component (time-domain)"""

        result = SuperpositionCurrent()
        for sub in self.sub.values():
            I = sub.get_I(name)
            result.add(I)
        result = result            
        return result

    def get_i(self, name):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def _get_Vd(self, Np, Nm=None):
        """This does not check nodes."""
        
        result = SuperpositionVoltage()
        for sub in self.sub.values():
            Vd = sub.get_Vd(Np, Nm)
            result.add(Vd)
        result = result.canonical()
        return result

    def get_Vd(self, Np, Nm=None):
        """Voltage drop between nodes (time-domain)"""

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
    
    

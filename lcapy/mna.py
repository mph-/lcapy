"""
This module implements modified nodal analysis (MNA).

Copyright 2014, 2015 Michael Hayes, UCECE
"""

from __future__ import division
from lcapy.core import cExpr, VV, II, s, sqrt
from lcapy.core import Matrix, Vector
import sympy as sym

# Note, all the maths is performed using sympy expressions and the
# values and converted to Expr when required.  This is more
# efficient and, more importantly, overcomes some of the wrapping
# problems which casues the is_real attribute to be dropped.

def namelist(elements):
    return ', '.join([elt for elt in elements])


class Mdict(dict):

    def __init__(self, branchdict, **assumptions):

        super(Mdict, self).__init__()
        self.branchdict = branchdict
        # Hack, should compute assumptions on fly.
        self.assumptions = assumptions

    def __getitem__(self, key):

        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key

        if key in self.branchdict:
            n1, n2 = self.branchdict[key]
            return VV((self[n1] - self[n2]).simplify(), **self.assumptions)

        return super(Mdict, self).__getitem__(key)

    def keys(self):

        return list(super(Mdict, self).keys()) + list(self.branchdict)


class MNA(object):
    """This class performs modified nodal analysis (MNA) on a netlist of
    components.  There are several variants:
    
    1. DC analysis if all the independent sources are DC.  The .V and .I
    methods return s-domain expressions with the dc assumption set.

    2. AC analysis if all the independent sources are AC.  The .V and .I
    methods return phasors.

    3. Initial value Laplace analysis if an L or C has an explicit
    initial value.  The .V and .I methods return s-domain expressions
    with no assumption sets; thus the time-domain results are only
    valid for t >= 0.

    4. General Laplace analysis.  If all the sources are causal and
    all the initial conditions are zero (explicitly or implicitly)
    then the time-domain results are causal.

    5. Mixed analysis.  This is not yet supported.  When the
    independent sources are mixed (say AC and DC), superpostion can be
    employed.

    """

    @property
    def is_causal(self):
        """Return True if all independent sources are causal"""

        independent_sources = self.independent_sources
        if independent_sources == {}:
            return not self.initial_value_problem

        for elt in self.independent_sources.values():
            if not elt.is_causal:
                return False
        return True

    @property
    def is_dc(self):
        """Return True if all independent sources are DC."""

        independent_sources = self.independent_sources
        if independent_sources == {}:
            return not self.initial_value_problem

        for elt in independent_sources.values():
            if not elt.is_dc:
                return False
        return True

    @property
    def is_ac(self):
        """Return True if all independent sources are AC."""

        independent_sources = self.independent_sources
        if independent_sources == {}:
            return not self.initial_value_problem

        for elt in independent_sources.values():
            if not elt.is_ac:
                return False
        return True

    @property
    def assumptions(self):

        assumptions = {}
        if self.is_ac:
            assumptions['ac'] = True
        if self.is_dc:
            assumptions['dc'] = True
        if self.is_causal:
            assumptions['causal'] = True
        return assumptions

    @property
    def zeroic(self):
        """Return True if the initial conditions for all components are zero"""

        for elt in self.elements.values():
            if not elt.zeroic:
                return False
        return True

    @property
    def initial_value_problem(self):
        """Return True if any components that allow initial conditions
        have them explicitly defined."""

        for elt in self.elements.values():
            if elt.hasic is None:
                continue
            if elt.hasic:
                return True

        return False

    @property
    def missing_ic(self):
        """Return components that allow initial conditions but do not have
        them explicitly defined"""

        return dict((key, elt) for key, elt in self.elements.items() if elt.hasic is False)

    @property
    def explicit_ic(self):
        """Return components that have explicitly defined initial conditions

        """

        return dict((key, elt) for key, elt in self.elements.items() if elt.hasic is True)

    @property
    def allow_ic(self):
        """Return components (L and C) that allow initial conditions"""

        return dict((key, elt) for key, elt in self.elements.items() if elt.hasic is not None)

    @property
    def noncausal_sources(self):
        """Return non-causal independent sources"""

        return dict((key, elt) for key, elt in self.elements.items() if elt.source and not elt.is_causal)

    @property
    def independent_sources(self):
        """Return independent sources (this does not include
        implicit sources due to initial conditions)"""

        return dict((key, elt) for key, elt in self.elements.items() if elt.source)

    @property
    def lnodes(self):
        """Determine linked nodes (both implicitly and explicitly
        connected)"""

        from copy import deepcopy

        # Start with implicitly linked nodes.
        lnodes = deepcopy(self.snodes)

        # Then augment with nodes connected by wires.
        for m, elt in enumerate(self.elements.values()):
            if elt.type not in ('W', ):
                continue

            n1, n2 = elt.nodes

            for key1, nodes in lnodes.items():
                if n1 in nodes:
                    break

            for key2, nodes in lnodes.items():
                if n2 in nodes:
                    break

            if key1 != key2:
                lnodes[key1].extend(lnodes.pop(key2))

        # Remove nodes that are not linked.
        pnodes = []
        for nodes in lnodes.values():
            if len(nodes) > 1:
                pnodes.append(nodes)

        return pnodes

    @property
    def node_map(self):
        """Determine mapping of nodes to unique nodes of the same potential"""

        if hasattr(self, '_node_map'):
            return self._node_map

        lnodes = self.lnodes

        node_map = {}
        for node in self.nodes:

            node_map[node] = node
            for nodes in lnodes:
                if node in nodes:
                    # Use first of the linked nodes unless '0' in list
                    if '0' in nodes:
                        node_map[node] = '0'
                    else:
                        # Sort nodes so 8 before 8_1 etc.
                        node_map[node] = sorted(nodes)[0]
                    break

        if '0' not in node_map:
            raise RuntimeError('Nothing connected to ground node 0')

        self._node_map = node_map
        return node_map

    @property
    def node_list(self):
        """Determine list of unique nodes"""

        if hasattr(self, '_node_list'):
            return self._node_list

        # Extract unique nodes.
        node_list = list(set(self.node_map.values()))
        # Ensure node '0' is first in the list.
        node_list.insert(0, node_list.pop(node_list.index('0')))

        self._node_list = node_list
        return node_list

    def _node_index(self, node):
        """Return node index; ground is -1"""
        return self.node_list.index(self.node_map[node]) - 1

    def _branch_index(self, cpt_name):

        try:
            index = self.unknown_branch_currents.index(cpt_name)
            return index
        except ValueError:
            raise ValueError('Unknown component name %s for branch current' % cpt_name)

    def _analyse(self):
        """Analyse network."""

        if hasattr(self, '_A'):
            return

        # Hack, to indirectly generate element list for network.
        if self.elements == {}:
            raise ValueError('No elements to analyse')

        # TODO: think this out.  When a circuit is converted
        # to a s-domain model we get Z (and perhaps Y) components.
        # We also loose the ability to determine the voltage
        # across a capacitor or inductor since they get split
        # into a Thevenin model and renamed.
        if hasattr(self, '_s_model'):
            raise RuntimeError('Cannot analyse s-domain model')
            
        if '0' not in self.node_map:
            raise RuntimeError('Nothing connected to ground node 0')

        # Determine if all components that allow initial conditions
        # have them explicitly defined.  In this case, we can only
        # provide solution for t >= 0.
        if self.is_causal and self.initial_value_problem and not self.zeroic:
            raise RuntimeError('Detected initial value problem that has causal sources!')

        if (not self.is_ac and not self.initial_value_problem
            and not self.is_causal and self.missing_ic != {}):
            print('Warning non-causal sources detected (%s)'
                  ' and initial conditions missing for %s;'
                  ' expect unexpected transient!' % (
                      namelist(self.noncausal_sources),
                      namelist(self.missing_ic)))

        # Determine which branch currents are needed.
        self.unknown_branch_currents = []

        for elt in self.elements.values():
            if elt.type in ('E', 'H', 'L', 'V'):
                self.unknown_branch_currents.append(elt.name)

        # Generate stamps.
        num_nodes = len(self.node_list) - 1
        num_branches = len(self.unknown_branch_currents)

        self._G = sym.zeros(num_nodes, num_nodes)
        self._B = sym.zeros(num_nodes, num_branches)
        self._C = sym.zeros(num_branches, num_nodes)
        self._D = sym.zeros(num_branches, num_branches)

        self._Is = sym.zeros(num_nodes, 1)
        self._Es = sym.zeros(num_branches, 1)

        assumptions = self.assumptions

        for elt in self.elements.values():
            elt.stamp(self, **assumptions)

        # Augment the admittance matrix to form A matrix.
        self._A = self._G.row_join(self._B).col_join(self._C.row_join(self._D))
        # Augment the known current vector with known voltage vector
        # to form Z vector.
        self._Z = self._Is.col_join(self._Es)

    def _solve(self):
        """Solve network."""

        if hasattr(self, '_V'):
            return
        self._analyse()

        # Solve for the nodal voltages
        try:
            Ainv = self._A.inv()
        except ValueError:
            raise ValueError(
                'The MNA A matrix is not invertible; some nodes may need'
                ' connecting with high value resistors, a voltage source'
                ' might be short-circuited, a current source might be'
                ' open-circuited.')

        results = sym.simplify(Ainv * self._Z)

        results = results.subs(self.context.symbols)
        assumptions = self.assumptions

        branchdict = {}
        for elt in self.elements.values():
            if elt.type == 'K':
                continue
            n1, n2 = self.node_map[elt.nodes[0]], self.node_map[elt.nodes[1]]
            branchdict[elt.name] = (n1, n2)

        self.context.switch()

        # Create dictionary of node voltages
        self._V = Mdict(branchdict, **assumptions)
        self._V['0'] = VV(0, **assumptions)
        for n in self.nodes:
            index = self._node_index(n)
            if index >= 0:
                self._V[n] = VV(results[index].simplify(), **assumptions)
            else:
                self._V[n] = VV(0, **assumptions)

        num_nodes = len(self.node_list) - 1

        # Create dictionary of branch currents through elements
        self._I = {}
        for m, key in enumerate(self.unknown_branch_currents):
            self._I[key] = II(results[m + num_nodes].simplify(), **assumptions)

        # Calculate the branch currents.  These should be lazily
        # evaluated as required.
        for elt in self.elements.values():
            if elt.type in ('R', 'C'):
                n1, n2 = self.node_map[
                    elt.nodes[0]], self.node_map[elt.nodes[1]]
                V1, V2 = self._V[n1], self._V[n2]
                I = ((V1 - V2 - elt.Vcpt) / elt.Zcpt)
                self._I[elt.name] = II(I.simplify(), **assumptions)

        self.context.restore()

    @property
    def A(self):
        """Return A matrix for MNA"""

        self._analyse()
        return Matrix(self._A)

    @property
    def Z(self):
        """Return Z vector for MNA"""

        self._analyse()
        return Vector(self._Z)

    @property
    def X(self):
        """Return X vector (of unknowns) for MNA"""

        self._analyse()

        V = ['V_' + node for node in self.node_list[1:]]
        I = ['I_' + branch for branch in self.unknown_branch_currents]
        return Vector(V + I)

    @property
    def V(self):
        """Return dictionary of s-domain node voltages indexed by node name"""

        self._solve()
        return self._V

    @property
    def I(self):
        """Return dictionary of s-domain branch currents indexed
        by component name"""

        self._solve()
        return self._I

    @property
    def Vd(self):
        """Return dictionary of s-domain branch voltage differences
        indexed by component name"""

        if hasattr(self, '_Vd'):
            return self._Vd

        self._solve()

        # This is a hack.  The assumptions should be recalculated
        # when performing operations on Expr types.
        assumptions = self.assumptions

        self._Vd = {}
        for elt in self.elements.values():
            if elt.type == 'K':
                continue
            n1, n2 = self.node_map[elt.nodes[0]], self.node_map[elt.nodes[1]]
            self._Vd[elt.name] = VV((self._V[n1] - self._V[n2]).simplify(),
                                    **assumptions)

        return self._Vd

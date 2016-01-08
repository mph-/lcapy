"""
This module implements modified nodal analysis (MNA).

Copyright 2014, 2015 Michael Hayes, UCECE
"""

from __future__ import division
from lcapy.core import cExpr, Vs, Is, s, sqrt
from lcapy.twoport import Matrix, Vector
import sympy as sym

# Note, all the maths is performed using sympy expressions and the
# values and converted to Expr when required.  This is more
# efficient and, more importantly, overcomes some of the wrapping
# problems which casues the is_real attribute to be dropped.

class Mdict(dict):

    def __init__(self, branchdir):

        super(Mdict, self).__init__()
        self.branchdir = branchdir

    def __getitem__(self, key):

        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key

        if key in self.branchdir:
            n1, n2 = self.branchdir[key]
            return self[n1] - self[n2]

        return super(Mdict, self).__getitem__(key)

    def keys(self):

        return super(Mdict, self).keys() + self.branchdir.keys()


class MNA(object):

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

            for key1, nodes in lnodes.iteritems():
                if n1 in nodes:
                    break

            for key2, nodes in lnodes.iteritems():
                if n2 in nodes:
                    break

            if key1 != key2:
                lnodes[key1].extend(lnodes.pop(key2))

        # Remove nodes that are not linked.
        pnodes = []
        for key, nodes in lnodes.iteritems():
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
            print('Nothing connected to ground node 0')
            node_map['0'] = '0'

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

        # TODO: think this out.  When a circuit is converted
        # to a s-domain model we get Z (and perhaps Y) components.
        # We also loose the ability to determine the voltage
        # across a capacitor or inductor since they get split
        # into a Thevenin model and renamed.
        if hasattr(self, '_s_model'):
            raise RuntimeError('Cannot analyse s-domain model')
            
        if '0' not in self.node_map:
            print('Nothing connected to ground node 0')
            self.nodes['0'] = None

        # Determine which branch currents are needed.
        self.unknown_branch_currents = []

        for key, elt in self.elements.iteritems():
            if elt.type in ('E', 'H', 'L', 'V'):
                self.unknown_branch_currents.append(key)

        # Generate stamps.
        num_nodes = len(self.node_list) - 1
        num_branches = len(self.unknown_branch_currents)

        self._G = sym.zeros(num_nodes, num_nodes)
        self._B = sym.zeros(num_nodes, num_branches)
        self._C = sym.zeros(num_branches, num_nodes)
        self._D = sym.zeros(num_branches, num_branches)

        self._Is = sym.zeros(num_nodes, 1)
        self._Es = sym.zeros(num_branches, 1)

        for elt in self.elements.values():
            elt.stamp(self)

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

        branchdir = {}
        for elt in self.elements.values():
            if elt.type == 'K':
                continue
            n1, n2 = self.node_map[elt.nodes[0]], self.node_map[elt.nodes[1]]
            branchdir[elt.name] = (n1, n2)

        # Create dictionary of node voltages
        self._V = Mdict(branchdir)
        self._V['0'] = Vs(0)
        for n in self.nodes:
            index = self._node_index(n)
            if index >= 0:
                self._V[n] = Vs(results[index])
            else:
                self._V[n] = Vs(0)

        num_nodes = len(self.node_list) - 1

        # Create dictionary of branch currents through elements
        self._I = {}
        for m, key in enumerate(self.unknown_branch_currents):
            self._I[key] = Is(results[m + num_nodes])

        # Calculate the branch currents.  These should be lazily
        # evaluated as required.
        for key, elt in self.elements.iteritems():
            if elt.type in ('R', 'C'):
                n1, n2 = self.node_map[
                    elt.nodes[0]], self.node_map[elt.nodes[1]]
                V1, V2 = self._V[n1], self._V[n2]
                I = ((V1 - V2 - elt.cpt.V) / elt.cpt.Z).simplify()
                self._I[elt.name] = Is(I)

    @property
    def A(self):
        """Return A matrix for MNA"""

        return Matrix(self._A)

    @property
    def Z(self):
        """Return Z vector for MNA"""

        return Vector(self._Z)

    @property
    def X(self):
        """Return X vector (of unknowns) for MNA"""

        V = ['V_' + node for node in self.node_list[1:]]
        I = ['I_' + branch for branch in self.unknown_branch_currents]
        return Vector(V + I)

    @property
    def V(self):
        """Return dictionary of s-domain node voltages indexed by node name"""

        if not hasattr(self, '_V'):
            self._solve()
        return self._V

    @property
    def I(self):
        """Return dictionary of s-domain branch currents indexed
        by component name"""

        if not hasattr(self, '_I'):
            self._solve()
        return self._I

    @property
    def Vd(self):
        """Return dictionary of s-domain branch voltage differences
        indexed by component name"""

        if hasattr(self, '_Vd'):
            return self._Vd

        self._Vd = {}
        for elt in self.elements.values():
            if elt.type == 'K':
                continue
            n1, n2 = self.node_map[elt.nodes[0]], self.node_map[elt.nodes[1]]
            self._Vd[elt.name] = Vs(sym.simplify(self.V[n1] - self.V[n2]))

        return self._Vd

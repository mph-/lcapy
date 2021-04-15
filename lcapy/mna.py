"""
This module implements modified nodal analysis (MNA).

Copyright 2014--2019 Michael Hayes, UCECE
"""

from __future__ import division
from .assumptions import Assumptions
from .vector import Vector
from .matrix import Matrix, matrix_inverse
from .sym import symsimplify
from .expr import ExprDict, expr
from .voltage import Vtype
from .current import Itype
from .systemequations import SystemEquations
import sympy as sym

# Note, all the maths is performed using sympy expressions and the
# values and converted to Expr when required.  This is more
# efficient and, more importantly, overcomes some of the wrapping
# problems which casues the is_real attribute to be dropped.

class Nodedict(ExprDict):

    def __getitem__(self, name):
        """Return node by name or number."""

        # If name is an integer, convert to a string.
        if isinstance(name, int):
            name = '%d' % name
        return super(Nodedict, self).__getitem__(name)


class Branchdict(ExprDict):
    pass
    

class MNAMixin(object):
    """This class performs modified nodal analysis (MNA) on a netlist of
    components.  There are several variants:
    
    1. DC analysis if all the independent sources are DC.  The .V and .I
    methods return DC expressions with the dc assumption set.

    2. AC analysis if all the independent sources are AC.  The .V and .I
    methods return phasors.

    3. Initial value Laplace analysis if an L or C has an explicit
    initial value.  The .V and .I methods return s-domain expressions
    with no assumption sets; thus the time-domain results are only
    valid for t >= 0.

    4. General Laplace analysis.  If all the sources are causal and
    all the initial conditions are zero (explicitly or implicitly)
    then the time-domain results are causal.

    5. Noise analysis.

    It is assumed that the use of this class uses superposition to
    solve problems with mixed independent sources, such as DC and AC.

    """

    def _invalidate(self):
        for attr in ('_A', '_Vdict', '_Idict'):
            if hasattr(self, attr):
                delattr(self, attr)

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
        # We also lose the ability to determine the voltage
        # across a capacitor or inductor since they get split
        # into a Thevenin model and renamed.
        if hasattr(self, '_s_model'):
            raise RuntimeError('Cannot analyse s-domain model')
            
        # Determine which branch currents are needed.
        self.unknown_branch_currents = []

        for elt in self.elements.values():
            if elt.need_branch_current:
                self.unknown_branch_currents.append(elt.name)
            if elt.need_extra_branch_current:
                self.unknown_branch_currents.append(elt.name + 'X')

        # Generate stamps.
        num_nodes = len(self.node_list) - 1
        num_branches = len(self.unknown_branch_currents)

        self._G = sym.zeros(num_nodes, num_nodes)
        self._B = sym.zeros(num_nodes, num_branches)
        self._C = sym.zeros(num_branches, num_nodes)
        self._D = sym.zeros(num_branches, num_branches)

        self._Is = sym.zeros(num_nodes, 1)
        self._Es = sym.zeros(num_branches, 1)

        # Iterate over circuit elements and fill in matrices.
        for elt in self.elements.values():
            elt._stamp(self)

        # Augment the admittance matrix to form A matrix.
        self._A = self._G.row_join(self._B).col_join(self._C.row_join(self._D))
        # Augment the known current vector with known voltage vector
        # to form Z vector.
        self._Z = self._Is.col_join(self._Es)

    def _solve(self):
        """Solve network."""
        
        if hasattr(self, '_Vdict'):
            return
        self._analyse()

        if '0' not in self.node_map:
            raise RuntimeError('Cannot solve: nothing connected to ground node 0')
        
        # Solve for the nodal voltages
        try:
            # The default method, Gaussian elimination, is the fastest
            # but hangs on some matrices with sympy-1.6.1
            # Comparative times for the testsuites are:
            # GE 66, ADJ 73, LU 76. 
            Ainv = matrix_inverse(self._A)
        except ValueError:
            comment = ''
            if self.kind == 'dc':
                comment = '  Check there is a DC path between all nodes.'
            raise ValueError(
"""The MNA A matrix is not invertible for %s analysis because:
1. there may be capacitors in series;
2. a voltage source might be short-circuited;
3. a current source might be open-circuited;
4. a dc current source is connected to a capacitor (use step current source).
5. part of the circuit is not referenced to ground
%s""" % (self.kind, comment))

        results = symsimplify(Ainv * self._Z)

        results = results.subs(self.context.symbols)

        branchdict = {}
        for elt in self.elements.values():
            if elt.type == 'K' or elt.ignore:
                continue
            n1, n2 = self.node_map[elt.nodenames[0]], self.node_map[elt.nodenames[1]]
            branchdict[elt.name] = (n1, n2)

        vtype = Vtype(self.kind)
        itype = Itype(self.kind)
        assumptions = Assumptions()
        if vtype.is_phasor_domain:
            assumptions.set('omega', self.kind)
        elif self.kind in ('s', 'ivp'):
            assumptions.set('ac', self.is_ac)
            assumptions.set('dc', self.is_dc)
            assumptions.set('causal', self.is_causal)            
        elif isinstance(self.kind, str) and self.kind[0] == 'n':
            assumptions.set('nid', self.kind)
       
        # Create dictionary of node voltages
        self._Vdict = Nodedict()
        self._Vdict['0'] = vtype(0, **assumptions)
        for n in self.nodes:
            index = self._node_index(n)
            if index >= 0:
                self._Vdict[n] = vtype(results[index], **assumptions).simplify()
            else:
                self._Vdict[n] = vtype(0, **assumptions)

        num_nodes = len(self.node_list) - 1

        # Create dictionary of branch currents through elements
        self._Idict = Branchdict()
        for m, key in enumerate(self.unknown_branch_currents):
            I = results[m + num_nodes]
            if key in self.elements and self.elements[key].is_source:
                I = -I
            self._Idict[key] = itype(I, **assumptions).simplify()

        # Calculate the branch currents.  These should be lazily
        # evaluated as required.
        for elt in self.elements.values():
            if elt.type in ('R', 'NR', 'C'):
                n1 = self.node_map[elt.nodenames[0]]
                n2 = self.node_map[elt.nodenames[1]]                
                V1, V2 = self._Vdict[n1], self._Vdict[n2]
                I = (V1.expr - V2.expr - elt.V0.expr) / elt.Z.expr
                self._Idict[elt.name] = itype(I, **assumptions).simplify()
            elif elt.type in ('I', ):
                self._Idict[elt.name] = elt.Isc

    @property
    def A(self):
        """Return A matrix for MNA"""

        self._analyse()
        return Matrix(self._A)

    @property
    def ZV(self):
        """Return Z vector for MNA"""

        self._analyse()
        return Vector(self._Z)

    @property
    def X(self):
        """Return X vector (of unknowns) for MNA"""

        self._analyse()
        
        V = [self.Vname('Vn%s' % node) for node in self.node_list[1:]]
        I = [self.Iname('I%s' % branch) for branch in self.unknown_branch_currents]
        return Vector(V + I)

    @property
    def Vdict(self):
        """Return dictionary of transform domain node voltages indexed by node
        name"""

        self._solve()
        return self._Vdict

    @property
    def Idict(self):
        """Return dictionary of transform domain branch currents indexed by
        component name"""

        self._solve()
        return self._Idict

    def matrix_equations(self, form='default', invert=False):
        """System of equations used to find the unknowns.

        Forms can be:
         A y = b
         b = A y
         Ainv b = y
         y = Ainv b

        If `invert` is True, evaluate the matrix inverse."""

        self._analyse()
        
        sys = SystemEquations(self._A, self._Z, self.X)        
        return sys.format(form, invert)        

    def equations(self, inverse=False):
        """System of equations used to find the unknowns.

        If inverse is True, evaluate the matrix inverse.

        This is for compatibility and is deprecated.  Use 
        matrix_equations instead."""

        return self.matrix_equations(invert=inverse)
    

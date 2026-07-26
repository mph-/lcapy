"""
This module performs loop analysis.  It is primarily for showing
the equations rather than evaluating them.

Copyright 2019--2026 Michael Hayes, UCECE

"""

# TODO, handle current sources and dependent sources.

from .circuitgraph import CircuitGraph
from .expr import expr, equation, ExprTuple
from .systemequations import SystemEquations
import sympy as sym


class LoopAnalysis(object):
    """This performs for loop analysis.  Currently, it uses mesh
    analysis and so is only applicable to circuits with a planar
    topology and without twoport components.

    Note, the loops are found using Networkx and since this uses
    unordered Python sets, different invocations can find different
    current loops for the same circuit.

    >>> from lcapy import Circuit, LoopAnalysis
    >>> cct = Circuit('''
    ... V1 1 0 {u(t)}; down
    ... R1 1 2; right=2
    ... L1 2 3; down=2
    ... W1 0 3; right
    ... W 1 5; up
    ... W 2 6; up
    ... C1 5 6; right=2
    ...''')

    To perform loop analysis in the time domain:

    >>> la = LoopAnalysis(cct)

    To display the equations found by applying KVL around each mesh:

    >>> la.mesh_equations().pprint()

    To display the system of equations (in matrix form) that needs to
    be solved:

    >>> la.matrix_equations().pprint()

    This only works for dc, ac, or laplace domains.  For example,

    >>> LoopAnalysis(cct.laplace()).matrix_equations().pprint()

    """

    # twoport components could be handled by splitting into
    # input and output oneports.  In the case of an ideal transformer,
    # an unknown voltage across the primary could be added to the mesh
    # equations.  The voltage across the secondary is the primary voltage
    # multiplied by N2 / N1.  An auxiliary equation is required
    # relating the mesh currents in the primary and secondary loops:
    # ib = -N1 / N2 ia.

    @classmethod
    def from_circuit(cls, cct):

        return cls(cct)

    def __init__(self, cct):

        self.cct = cct
        self.cg = CircuitGraph.from_circuit(cct)

        self.kind = cct.kind
        if cct.kind == 'super':
            self.kind = self.cct._analysis_kind()

        self._equations = self._make_equations()

        self._unknowns = self.mesh_currents()
        self._y = matrix(self._unknowns)

    def loops(self):
        """Return list of basis loops.  Note, the loops can vary for different
        invocations of the LoopAnalysis class.

        See also loops_by_cpt() and loops_by_cpt_name().
        """

        return self.cg.circuit_loops()

    def loops_by_cpt(self):
        """Return list of loops specified by cpt.

        See also loops() and loops_by_cpt_name().
        """

        return self.cg.basis_loops_by_cpt_name()

    def loops_by_cpt_name(self):
        """Return list of loops specified by cpt name.

        See also loops() and loops_by_cpt().
        """

        return self.cg.basis_loops_by_cpt_name()

    def loops_with_cpt_name(self, name):
        """Return list of loops containing cpt name.
        """

        loops = []
        for loop in self.loops():
            if loop.has_cpt_name(name):
                loops.append(loop)
                continue
        return loops

    def edges_with_cpt_name(self, name):
        """Return list of edges containing cpt name.
        """

        edges = []
        for loop in self.loops():
            for edge in loop:
                if edge.has_cpt_name(name):
                    edges.append(edge)
                    continue
        return edges

    @property
    def num_loops(self):

        return len(self.cg.basis_loops)

    def mesh_currents(self):
        """Return list of mesh current names."""

        if not self.cg.is_planar:
            raise ValueError('Circuit topology is not planar')

        loops = self.loops()
        Nloops = len(loops)

        mesh_currents = ExprList(
            [Iname('I_%d' % (m + 1), self.kind) for m in range(Nloops)])
        return mesh_currents

    def _add_mesh_currents(self, edge, loop, loops, mesh_currents):
        """Sum the currents flowing through the specified edge for all
        loops."""

        net_current = Itype(self.kind)(0)

        cpt_name = edge.cpt.name

        for n, loop2 in enumerate(loops):

            if not loop2.has_cpt_name(cpt_name):
                continue
            edge2 = loop2.edge_by_cpt_name(cpt_name)

            if edge.node1 == edge2.node1 and edge.node2 == edge2.node2:
                polarity = 1
            elif edge.node1 == edge2.node2 and edge.node2 == edge2.node1:
                polarity = -1
            else:
                raise ValueError('Unmatched edge')

            net_current += polarity * mesh_currents[n]

        return net_current

    def _process_loop(self, loop, mesh_current, loops, mesh_currents):
        """Generate mesh equation for specified loop."""

        result = Vtype(self.kind)(0)

        for edge in loop:
            cpt = edge.cpt
            # Skip dummy wires.
            if cpt is None:
                continue

            if cpt.is_current_source:
                raise ValueError('TODO: handle current source in loop')
            elif cpt.is_dependent_source:
                raise ValueError('Dependent sources not handled yet')

            if cpt.is_voltage_source:
                v = cpt.cpt.voltage_equation(mesh_current, self.kind)
                if edge.is_flipped:
                    v = -v
            else:
                current = self._add_mesh_currents(edge, loop, loops,
                                                  mesh_currents)
                v = -cpt.cpt.voltage_equation(current, self.kind)

            result += v

        return result

    def _make_equations(self):

        if not self.cg.is_planar:
            raise ValueError('Circuit topology is not planar')

        loops = self.loops()
        Nloops = len(loops)

        mesh_currents = [Iname('I_%d' % (m + 1), self.kind)
                         for m in range(Nloops)]
        equations = {}

        # For each loop, generate the mesh equation.
        for m, loop in enumerate(loops):

            result = self._process_loop(loop, mesh_currents[m],
                                        loops, mesh_currents)
            equations[mesh_currents[m]] = (result, expr(0))

        return equations

    def mesh_equations_list(self):
        """Return mesh equations as a list."""

        result = ExprList()

        for current, (lhs, rhs) in self._equations.items():
            result.append(equation(lhs, rhs))
        return result

    def mesh_equations(self):
        """Return mesh equations as a dict keyed by the mesh current."""

        result = ExprDict()

        for current, (lhs, rhs) in self._equations.items():
            result[current] = equation(lhs, rhs)
        return result

    def _analyse(self):

        if self.kind in ('t', 'time'):
            raise ValueError(
                'Cannot put time domain equations into matrix form.  '
                'Convert to dc, ac, or laplace domain first.')

        subsdict = {}
        for m, i in enumerate(self._unknowns):
            subsdict[i.expr] = 'X_X%d' % m

        exprs = []
        for node, (lhs, rhs) in self._equations.items():
            lhs = lhs.subs(subsdict).expr.expand()
            rhs = rhs.subs(subsdict).expr.expand()
            exprs.append(lhs - rhs)

        y = []
        for y1 in self._y:
            y.append(y1.subs(subsdict).expr)

        A, b = sym.linear_eq_to_matrix(exprs, *y)

        y = [y1.expr for y1 in self._y]
        return SystemEquations(A, b, y)

    @property
    def A(self):
        """Return A matrix where A y = b."""

        if not hasattr(self, '_sys'):
            self._sys = self._analyse()
        return matrix(self._sys.A)

    @property
    def b(self):
        """Return b vector where A y = b."""

        if not hasattr(self, '_sys'):
            self._sys = self._analyse()
        return matrix(self._sys.b)

    @property
    def y(self):
        """Return y vector where A y = b."""
        return self._y

    def matrix_equations(self, form='default', invert=False):
        """Return the equations in matrix form.

        Forms can be:
         'default'
         'A y = b'
         'b = A y'
         'Ainv b = y'
         'y = Ainv b'

        If `invert` is True, evaluate the matrix inverse."""

        if not hasattr(self, '_sys'):
            self._sys = self._analyse()
        return self._sys.format(form, invert)

    @property
    def unknowns(self):
        """Return tuple of the unknown voltages"""

        return ExprTuple(self.y)

    def solve_laplace(self):
        """Determine the unknown voltages using Laplace transforms and
        return as a dict"""

        from .sexpr import s

        unknowns = self.unknowns(s)
        return self.mesh_equations()(s).solve(unknowns)

    def pdb(self):
        """Enter the python debugger."""

        import pdb
        pdb.set_trace()
        return self


from .expr import ExprList, ExprDict, expr  # nopep8
from .current import Iname, Itype  # nopep8
from .voltage import Vtype  # nopep8
from .matrix import matrix  # nopep8

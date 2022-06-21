"""
This module performs loop analysis.  It is primarily for showing
the equations rather than evaluating them.

Copyright 2019--2022 Michael Hayes, UCECE

"""

# TODO, handle current sources and dependent sources.

from .circuitgraph import CircuitGraph
from .expr import expr, equation, ExprTuple
from .systemequations import SystemEquations
import sympy as sym


class LoopAnalysis(object):
    """This is an experimental class for loop analysis.  Currently,
    it uses mesh analysis and so is only applicable to circuits with
    a planar topology.

    The API is likely to change since different invocations find
    different current loops.

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

    """

    @classmethod
    def from_circuit(cls, cct):

        return cls(cct)

    def __init__(self, cct):

        self.cct = cct
        self.cg = CircuitGraph.from_circuit(cct)

        self.kind = self.cct.kind
        if self.kind == 'super':
            self.kind = 'time'

        self._equations = self._make_equations()

        self._unknowns = self.mesh_currents()
        self._y = matrix(self._unknowns)

    def loops(self):
        """Return list of loops.  Note, the loops can vary for different
        invocations of the LoopAnalysis class."""

        return self.cg.loops()

    def loops_by_cpt_name(self):
        """Return list of loops specified by cpt name."""

        return self.cg.loops_by_cpt_name()

    @property
    def num_loops(self):

        return len(self.cg.loops)

    def mesh_currents(self):

        if not self.cg.is_planar:
            raise ValueError('Circuit topology is not planar')

        loops = self.loops()
        Nloops = len(loops)

        mesh_currents = ExprList(
            [Iname('I_%d' % (m + 1), self.kind) for m in range(Nloops)])
        return mesh_currents

    def _make_equations(self):

        if not self.cg.is_planar:
            raise ValueError('Circuit topology is not planar')

        if self.cct.is_superposition:
            raise ValueError(
                'Circuit has a mixture of ac/dc/transient sources')

        # This is work in progress to do mesh analysis.

        loops = self.loops()
        Nloops = len(loops)

        mesh_currents = [Iname('I_%d' % (m + 1), self.kind)
                         for m in range(Nloops)]
        equations = {}

        for m, loop in enumerate(loops):

            result = Vtype(self.kind)(0)

            loop1 = loop.copy()
            loop1.append(loop1[0])
            for j in range(len(loop1) - 1):

                elt = self.cg.component(loop1[j], loop1[j + 1])
                if elt is None:
                    continue

                if elt.is_current_source:
                    raise ValueError('TODO: handle current source in loop')

                # Map node names to equipotential node names.
                nodenames = [self.cct.node_map[nodename]
                             for nodename in elt.nodenames]

                is_reversed = nodenames[0] == loop1[j] and nodenames[1] == loop1[j + 1]

                if is_reversed:
                    current = -mesh_currents[m]
                else:
                    current = mesh_currents[m]

                for n, loop2 in enumerate(loops):

                    if loop2 == loop:
                        continue

                    loop2 = loop2.copy()
                    loop2.append(loop2[0])

                    if nodenames[0] not in loop2 or nodenames[1] not in loop2:
                        continue

                    # Determine polarity of cpt.
                    for l in range(len(loop2) - 1):
                        if (nodenames[0] == loop2[l] and
                                nodenames[1] == loop2[l + 1]):
                            current -= mesh_currents[n]
                            break
                        elif (nodenames[1] == loop2[l] and
                              nodenames[0] == loop2[l + 1]):
                            current += mesh_currents[n]
                            break

                v = elt.cpt.voltage_equation(current, self.kind)
                if elt.is_voltage_source and is_reversed:
                    v = -v
                result += v

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

        if self.kind == 'time':
            raise ValueError(
                'Cannot put time domain equations into matrix form')

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


from .expr import ExprList, ExprDict, expr  # nopep8
from .current import Iname  # nopep8
from .voltage import Vtype  # nopep8
from .matrix import matrix  # nopep8

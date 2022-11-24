"""
This module performs state-space analysis.

Copyright 2021--2022 Michael Hayes, UCECE

"""

from .mnacpts import I, V, C, L
from .expr import expr
from .matrix import Matrix
from .tmatrix import TimeDomainMatrix
from .voltage import voltage
from .current import current
from .statespace import StateSpace
from .circuitgraph import CircuitGraph
from .sym import sympify
import sympy as sym
from warnings import warn

__all__ = ('StateSpaceMaker', )

# TODO
# 1. Use a better Matrix class that preserves the class of each
# element, where possible.
# 2. Handle namespaces

# Independent sources
# V1 1 0
# V1 1 0 {10 * u(t)}
# V1 1 0 ac 10
#
# In each case, we use v1(t) as the independent source when determining
# the A, B, C, D matrices.  We can then substitute the known value at
# the end.

# There is a generalised (but less common) state-space representation:
#
# dx/dt = A x + B u + C du/dt
# y = D x + E u + F du/dt
#
# This is required when have a current source in series with an inductor;
# here the inductor current is not a state variable since it is equivalent
# to an input.  Moreover, the output equation for the inductor voltage
# requires the derivative of the input current.


def _hack_vars(exprs):
    """Substitute i_Canon1(t) with i_C(t) etc. provided
    there is no i_Canon2(t)."""

    for m, expr1 in enumerate(exprs):
        for c in ('i_V', 'i_C', 'i_L', 'v_C'):
            sym1 = sympify(c + 'anon1(t)')
            sym2 = sympify(c + 'anon2(t)')
            if expr1.has(sym1) and not expr1.has(sym2):
                expr1 = expr1.subs(sym1, sympify(c + '(t)'))
                exprs[m] = expr1


class StateSpaceMaker(object):
    """This converts a circuit to state-space representation.

    It can take a long time for a symbolic circuit with many reactive
    components.

    The currents through inductors and the voltage across capacitors
    are chosen as the state variables.  However, this fails when a
    current source is in series with an inductor.  In this case, the
    inductor current is not a state variable since it is equivalent to
    an input.  Moreover, the output equation for the inductor voltage
    requires the derivative of the input current.

    This does not (yet) look for degenerate circuits.  These are
    circuits with a loop consisting only of voltage sources and/or
    capacitors, or a cut set consisting only of current sources and/or
    inductors.  One hack is to call simplify() first to remove series
    inductors and parallel capacitors.

    """

    @classmethod
    def from_circuit(cls, cct, node_voltages=None, branch_currents=None):
        """`node_voltages` is a list of node names to use as voltage outputs.
        If `None` use all the unique node names.  Use `()`
        if want no branch currents.

        `branch_currents` is a list of component names to use as
        current outputs.  If `None` use all the components.  Use `()`
        if want no branch currents."""

        if node_voltages is None:
            node_voltages = cct.node_list

        if branch_currents is None:
            branch_currents = cct.branch_list

        if node_voltages == [] and branch_currents == []:
            raise ValueError('State-space: no outputs')

        inductors = []
        capacitors = []
        independent_current_sources = []
        independent_voltage_sources = []

        # Determine state variables (current through inductors and
        # voltage across acapacitors) and replace inductors with
        # current sources and capacitors with voltage sources.
        sscct = cct._new()
        cpt_map = {}

        for key, elt in cct.elements.items():
            ssnet = elt._ss_model()
            sselt = sscct._add(ssnet)
            name = elt.name
            cpt_map[name] = sselt.name

            if elt.is_inductor:
                if sselt.name in cct.elements:
                    raise ValueError(
                        'Name conflict %s, either rename the component or improve the code!' % sselt.name)

                inductors.append(elt)
            elif elt.is_capacitor:
                if sselt.name in cct.elements:
                    raise ValueError(
                        'Name conflict %s, either rename the component or improve the code!' % sselt.name)
                capacitors.append(elt)
            elif elt.is_independent_current_source:
                independent_current_sources.append(elt)
            elif elt.is_independent_voltage_source:
                independent_voltage_sources.append(elt)

        independent_sources = independent_voltage_sources + independent_current_sources

        # Build circuit graph and check if have inductor in series with
        # current source.
        if independent_current_sources != [] and inductors != []:
            cg = CircuitGraph(cct)
            for elt in independent_current_sources:
                for name in cg.in_series(elt.name):
                    if cct[name].is_inductor:
                        raise ValueError(
                            'Cannot create state-space model since have inductor %s in series with independent current source %s' % (name, elt.name))

        cct = cct
        sscct = sscct
        # sscct can be analysed in the time domain since it has no
        # reactive components.  However, for large circuits
        # this can take a long time due to inversion of the MNA matrix.

        try:
            # Analyse node voltages and branch currents.
            sscct[0].V
        except ValueError as e:
            reasons = []
            if len(inductors) > 0:
                reasons.append(
                    'Check for cut set consisting only of current sources and/or inductors.')
            if len(capacitors) > 0:
                reasons.append(
                    'Check for a loop consisting of voltage sources and/or capacitors.')

            raise ValueError("State-space analysis failed.\n%s\n     %s" % (
                e, '\n    '.join(reasons)))

        dotx_exprs = []
        statevars = []
        statenames = []
        initialvalues = []
        for elt in inductors + capacitors:
            name = cpt_map[elt.name]

            if isinstance(elt, L):
                # Inductors  v = L di/dt  so need v across the L
                expr = sscct[name].v / elt.cpt.L
                var = -sscct[name].isc
                x0 = elt.cpt.i0
            else:
                # Capacitors  i = C dv/dt  so need i through the C
                # The current is negated since it is from a source V_Cx
                expr = -sscct[name].i / elt.cpt.C
                var = sscct[name].voc
                x0 = elt.cpt.v0

            dotx_exprs.append(expr)
            statevars.append(var)
            statenames.append(name)
            initialvalues.append(x0)

        statesyms = sympify(statenames)

        # Determine independent sources.
        sources = []
        sourcevars = []
        sourcenames = []
        for elt in independent_sources:
            name = cpt_map[elt.name]

            if isinstance(elt, V):
                expr = elt.cpt.voc
                var = sscct[name].voc
            else:
                expr = elt.cpt.isc
                var = sscct[name].isc

            sources.append(expr)
            sourcevars.append(var)
            sourcenames.append(name)

        sourcesyms = sympify(sourcenames)

        # linear_eq_to_matrix expects only Symbols and not AppliedUndefs,
        # so substitute.
        subsdict = {}
        for var, sym1 in zip(statevars, statesyms):
            subsdict[var.expr] = sym1
        for var, sym1 in zip(sourcevars, sourcesyms):
            subsdict[var.expr] = sym1

        for m, expr in enumerate(dotx_exprs):
            dotx_exprs[m] = expr.subs(subsdict).expr.expand()

        if dotx_exprs == []:
            warn('State-space: no state variables found')
        if sourcesyms == []:
            warn('State-space: no independent sources found')

        if dotx_exprs != []:
            A, b = sym.linear_eq_to_matrix(dotx_exprs, *statesyms)
        else:
            A = sym.zeros(len(dotx_exprs), len(dotx_exprs))

        if sourcesyms != [] and dotx_exprs != []:
            B, b = sym.linear_eq_to_matrix(dotx_exprs, *sourcesyms)
        else:
            B = sym.zeros(len(dotx_exprs), len(sources))

        # Determine output variables.
        yexprs = []
        y = []

        for node in node_voltages:
            if node == '0':
                continue
            yexprs.append(sscct[node].v.subs(subsdict).expand().expr)
            # Note, this can introduce a name conflict
            y.append(voltage('v_%s(t)' % node))

        for name in branch_currents:
            # Perhaps ignore L since the current through it is a
            # state variable?
            name2 = cpt_map[name]
            yexprs.append(sscct[name2].i.subs(subsdict).expand().expr)
            y.append(current('i_%s(t)' % name))

        if statesyms != [] and yexprs != []:
            C, b = sym.linear_eq_to_matrix(yexprs, *statesyms)
        else:
            C = sym.zeros(len(yexprs), len(statesyms))

        if sourcesyms != [] and yexprs != []:
            D, b = sym.linear_eq_to_matrix(yexprs, *sourcesyms)
        else:
            D = sym.zeros(len(yexprs), len(sourcesyms))

        # Rewrite vCanon1(t) as vC(t) etc if appropriate.
        _hack_vars(statevars)
        _hack_vars(sources)

        # Note, Matrix strips the class from each element...
        x = TimeDomainMatrix(statevars)

        x0 = Matrix(initialvalues)

        u = TimeDomainMatrix(sources)

        # Perhaps could use v_R1(t) etc. as the output voltages?
        y = TimeDomainMatrix(y)

        A = Matrix(A)
        B = Matrix(B)
        C = Matrix(C)
        D = Matrix(D)

        return StateSpace(A, B, C, D, u, y, x, x0)

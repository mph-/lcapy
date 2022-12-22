"""
This module implements an experimental time-stepping simulation.

Copyright 2021--2022 Michael Hayes, UCECE
"""

from .sym import tsym, symbol_map
from .symbols import oo
from warnings import warn

__all__ = ('Simulator', )

# The companion circuits use a Thevenin model.  This simplifies
# determination of the current through reactive components.  If the
# currents are not required, the Norton model would be faster since
# fewer nodes are needed and so the matrices are smaller.

# TODO:
# 1. handle initial values
# 2. offset correction


class SimulatedComponent(object):

    def __init__(self, cpt, v1_index, v2_index, v3_index, i_index):

        self.nodes = cpt.nodenames
        self.name = cpt.name
        self.Reqname = 'R%seq' % cpt.name
        self.Veqname = 'V%seq' % cpt.name
        self.Reqsym = symbol_map(self.Reqname)
        self.Veqsym = symbol_map(self.Veqname)
        self.v1_index = v1_index
        self.v2_index = v2_index
        # This is the dummy node required for the Thevenin companion circuit.
        self.v3_index = v3_index
        self.i_index = i_index

    def subsdict(self, n, dt, v1, v2, i):
        """Create a dictionary of substitutions."""

        geq = self.geq(n, dt, v1, v2, i)
        veq = self.veq(n, dt, v1, v2, i)

        return {self.Reqsym: 1 / geq, self.Veqsym: veq}

    def stamp(self, A, Z, num_nodes, n, dt, v1, v2, i):

        geq = self.geq(n, dt, v1, v2, i)
        veq = self.veq(n, dt, v1, v2, i)

        n1, n2 = self.v1_index, self.v3_index

        if n1 >= 0 and n2 >= 0:
            A[n1, n2] -= geq
            A[n2, n1] -= geq
        if n1 >= 0:
            A[n1, n1] += geq
        if n2 >= 0:
            A[n2, n2] += geq

        m = self.i_index + num_nodes
        Z[m] += veq


class SimulatedCapacitor(SimulatedComponent):

    def __init__(self, C, v1_index, v2_index, v3_index, i_index):

        super(SimulatedCapacitor, self).__init__(C, v1_index, v2_index,
                                                 v3_index, i_index)
        self.Cval = C.C.expr


class SimulatedInductor(SimulatedComponent):

    def __init__(self, L, v1_index, v2_index, v3_index, i_index):

        super(SimulatedInductor, self).__init__(L, v1_index, v2_index,
                                                v3_index, i_index)
        self.Lval = L.L.expr


class SimulatedCapacitorTrapezoid(SimulatedCapacitor):

    def geq(self, n, dt, v1, v2, i):

        return (2 * self.Cval) / dt

    def veq(self, n, dt, v1, v2, i):

        if n < 1:
            return 0

        v = v1[n - 1] - v2[n - 1]

        geq = (2 * self.Cval) / dt
        veq = v + i[n - 1] / geq
        return veq


class SimulatedInductorTrapezoid(SimulatedInductor):

    def geq(self, n, dt, v1, v2, i):

        return dt / (2 * self.Lval)

    def veq(self, n, dt, v1, v2, i):

        if n < 1:
            return 0

        v = v1[n - 1] - v2[n - 1]

        geq = dt / (2 * self.Lval)
        veq = -v - i[n - 1] / geq
        return veq


class SimulatedCapacitorBackwardEuler(SimulatedCapacitor):

    def geq(self, n, dt, v1, v2, i):

        return self.Cval / dt

    def veq(self, n, dt, v1, v2, i):

        if n < 1:
            return 0

        return v1[n - 1] - v2[n - 1]


class SimulatedInductorBackwardEuler(SimulatedInductor):

    def geq(self, n, dt, v1, v2, i):

        return dt / self.Lval

    def veq(self, n, dt, v1, v2, i):

        geq = dt / self.Lval
        veq = -i[n - 1] / geq
        return veq


class SimulationResultsNode(object):

    def __init__(self, v):
        self.v = v


class SimulationResultsCpt(object):

    def __init__(self, v, i):
        self.v = v
        self.i = i


class SimulationResults(object):

    def __init__(self, tv, cct, r_model, node_list, branch_list):

        from numpy import zeros

        self.t = tv
        self.cct = cct
        self.r_model = r_model

        N = len(tv)

        # MNA calculates the node voltages plus currents through
        # inductors and voltage sources.  The currents through the
        # other elements can be found from the voltage difference
        # divided by the element resistance.

        self.num_nodes = len(node_list) - 1
        self.num_branches = len(branch_list)

        self.node_voltages = zeros((self.num_nodes + 1, N))
        self.branch_currents = zeros((self.num_branches, N))

    def __getitem__(self, name):
        """Return element or node by name."""

        cct = self.cct

        # If name is an integer, convert to a string.
        if isinstance(name, int):
            name = '%d' % name

        if name in cct.nodes:
            return SimulationResultsNode(self._node_voltage_get(name))

        if name in cct._elements:
            return SimulationResultsCpt(self._cpt_voltage_get(name),
                                        self._cpt_current_get(name))

        raise AttributeError('Unknown element or node name %s' % name)

    def __getattr__(self, attr):
        """Return element or node by name.
        This gets called if there is no explicit attribute attr for
        this instance.  This is primarily for accessing elements and
        non-numerical node names.  It also gets called if the called
        attr throws an AttributeError exception.  The annoying thing
        is that hasattr uses getattr and checks for an exception.

        """

        return self.__getitem__(attr)

    def _node_voltage_get(self, n):

        index = self.r_model.mna._node_index(n)
        # NB, node_voltages is zero for index = -1
        return self.node_voltages[index]

    def _cpt_voltage_get(self, cptname):

        cpt = self.cct.elements[cptname]

        v1 = self._node_voltage_get(cpt.nodenames[0])
        v2 = self._node_voltage_get(cpt.nodenames[1])
        return v1 - v2

    def _cpt_current_get(self, cptname):

        try:
            index = self.r_model.mna._branch_index(cptname)
            return self.branch_currents[index]
        except:
            cpt = self.cct._elements[cptname]
            if cpt.is_capacitor or cpt.is_inductor:
                # For a capacitor we can find the current through the
                # companion resistor or voltage source.
                return self._cpt_current_get('V%seq' % cptname)

            Vd = self._cpt_voltage_get(cptname)
            if cpt.is_resistor or cpt.is_conductor:
                return Vd / float(cpt.R.expr)

            # Need to determine resistance of the cpt
            raise ValueError('FIXME')

    @property
    def V(self, node):
        """Node voltage with respect to ground."""

        return self._node_voltage_get(node)


class Simulator(object):

    def __init__(self, cct):
        """Create simulation object for the circuit specified by `cct`.

        All the symbolic circuit component values need to be replaced
        with numerical values (using the subs method) except for
        functions of t, such as Heaviside(t).

        Here's an example of use:

        cct = Circuit('circuit.sch')
        sim = Simulator(cct)
        t = np.linspace(0, 1, 100)
        results = sim(t)

        plot(t, results.C1.v)

        """

        self.cct = cct

        # Companion resistor model
        self.r_model = cct.r_model().subcircuits['time']

    def _step(self, foo, n, tv, results):

        from numpy import array, linalg, dot

        # Substitute values into the MNA A matrix and Z vector,
        # then perform numerical inversion of the A matrix.
        # Alternatively, a faster way, if the A matrix was small,
        # would be to symbolically invert the A matrix and
        # then substitute values.

        if n == 0:
            if not self.cct.is_IVP:
                # Initial voltages and currents all zero.
                return
            p_model = self.cct.pre_initial_model()
            # Evaluate model and copy node voltages and branch currents...

            return

        dt = tv[n] - tv[n - 1]
        subsdict = {tsym: tv[n]}

        Zsym = self.Zsym
        Zsym = Zsym.subs(subsdict)
        Z = array(Zsym).astype(float).squeeze()

        if n == 1 and Zsym.free_symbols != set():
            raise ValueError(
                'Undefined symbols %s in Z vector; use subs to replace with numerical values' % Zsym.free_symbols)

        # Ensure have a copy.
        A = self.A + 0

        for cpt in self.reactive_cpts:

            # NB, node_voltages is zero for index = -1
            v1 = results.node_voltages[cpt.v1_index]
            v2 = results.node_voltages[cpt.v2_index]
            i = results.branch_currents[cpt.i_index]

            cpt.stamp(A, Z, results.num_nodes, n, dt, v1, v2, i)

        Ainv = linalg.inv(A)

        results1 = dot(Ainv, Z)

        num_nodes = results.num_nodes
        results.node_voltages[0:num_nodes, n] = results1[0:num_nodes]
        results.branch_currents[:, n] = results1[num_nodes:]

    def __call__(self, tv, integrator='trapezoid'):
        """Numerically evaluate circuit using time-stepping numerical
        integration at the vector of times specified by `tv`.

        Currently the only supported integration methods are
        trapezoidal and backward-euler (others would be trivial to
        add).  The trapezoidal integration method is the default since
        it is accurate but it can be unstable producing some
        oscillations.  Unfortunately, there is no ideal numerical
        integration method and there is always a tradeoff between
        accuracy and stability.

        """

        from numpy import array

        if integrator == 'trapezoid':
            Ccls = SimulatedCapacitorTrapezoid
            Lcls = SimulatedInductorTrapezoid
        elif integrator == 'backward-euler':
            Ccls = SimulatedCapacitorBackwardEuler
            Lcls = SimulatedInductorBackwardEuler
        else:
            raise ValueError('Unknown integrator ' + integrator)

        r_model = self.r_model

        Asubsdict = {}
        Zsubsdict = {}
        self.reactive_cpts = []
        for key, elt in self.cct.elements.items():
            if not (elt.is_inductor or elt.is_capacitor):
                continue

            if not elt.has_ic:
                warn('Initial conditions for %s ignored' % elt.name)

            v1_index = r_model.mna._node_index(elt.nodenames[0])
            v2_index = r_model.mna._node_index(elt.nodenames[1])
            i_index = r_model.mna._branch_index('V%seq' % elt.name)
            relt = self.r_model.elements['R%seq' % elt.name]
            v3_index = r_model.mna._node_index(relt.nodenames[1])

            if elt.is_inductor:
                cls = Lcls
            else:
                cls = Ccls

            simcpt = cls(elt, v1_index, v2_index, v3_index, i_index)
            self.reactive_cpts.append(simcpt)

            Asubsdict[simcpt.Reqsym] = oo
            Zsubsdict[simcpt.Veqsym] = 0

        # Remove 1 / Req entries
        Asym = r_model.mna._A.subs(Asubsdict)
        # Remove Veq entries
        Zsym = r_model.mna._Z.subs(Zsubsdict)

        self.Asym = Asym
        self.Zsym = Zsym

        if Asym.free_symbols != set():
            raise ValueError(
                'Undefined symbols %s in A matrix; use subs to replace with numerical values' % Asym.free_symbols)

        # Convert to numpy ndarray
        self.A = array(Asym).astype(float)

        results = SimulationResults(tv, self.cct, r_model,
                                    r_model.node_list,
                                    r_model.mna.unknown_branch_currents)

        for n, t1 in enumerate(tv):
            self._step(r_model, n, tv, results)

        return results

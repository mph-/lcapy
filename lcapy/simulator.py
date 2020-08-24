"""
This module implements an experimental time-stepping simulation.

Copyright 2020 Michael Hayes, UCECE
"""

from numpy import zeros, array, float, linalg, dot
from .sym import tsym, symbol_map

__all__ = ('Simulator', )

# The companion circuits use a Thevenin model.  This simplifies
# determination of the current through reactive components.  If the
# currents are not required, the Norton model would be faster since
# fewer nodes are needed and so the matrices are smaller.

# TODO:  speed up subs
# Replace 1 / R_eq with g_eq in A matrix.
# Record row, col, sign for each g_eq
# Replace g_eq with 0
# When time-stepping, add g_eq values into A matrix using stored info.

class SimulatedComponent(object):

    def __init__(self, cpt, v1_index, v2_index, i_index):
        
        self.nodes = cpt.nodes
        self.name = cpt.name
        self.Reqname = 'R%seq' % cpt.name
        self.Veqname = 'V%seq' % cpt.name
        self.Reqsym = symbol_map(self.Reqname)
        self.Veqsym = symbol_map(self.Veqname)                
        self.v1_index = v1_index
        self.v2_index = v2_index
        self.i_index = i_index
        

class SimulatedCapacitor(SimulatedComponent):

    def __init__(self, C, v1_index, v2_index, i_index):

        super (SimulatedCapacitor, self).__init__(C, v1_index, v2_index,
                                                  i_index)
        self.Cval = C.C.expr

    
class SimulatedInductor(SimulatedComponent):

    def __init__(self, L, v1_index, v2_index, i_index):

        super (SimulatedInductor, self).__init__(L, v1_index, v2_index,
                                                 i_index)
        self.Lval = L.L.expr

        
class SimulatedCapacitorTrapezoid(SimulatedCapacitor):

    def subsdict(self, n, dt, v1, v2, i):
        """Create a dictionary of substitutions."""

        v = v1[n - 1] - v2[n - 1]

        C = self.Cval

        Req = dt / (2 * C)
        veq = v + i[n - 1] * dt / (2 * C)

        return {self.Reqsym:Req, self.Veqsym:veq}    


class SimulatedInductorTrapezoid(SimulatedInductor):

    def subsdict(self, n, dt, v1, v2, i):
        """Create a dictionary of substitutions."""        

        v = v1[n - 1] - v2[n - 1]
        
        L = self.Lval
        
        Req = 2 * L / dt
        veq = -v - 2 * L * i[n - 1] / dt
        return {self.Reqsym:Req, self.Veqsym:veq}


class SimulatedCapacitorBackwardEuler(SimulatedCapacitor):

    def subsdict(self, n, dt, v1, v2, i):
        """Create a dictionary of substitutions."""

        v = v1[n - 1] - v2[n - 1]
        
        C = self.Cval     

        Req = dt / C
        veq = v

        return {self.Reqsym:Req, self.Veqsym:veq}    


class SimulatedInductorBackwardEuler(SimulatedInductor):

    def subsdict(self, n, dt, v1, v2, i):
        """Create a dictionary of substitutions."""        

        v = v1[n - 1] - v2[n - 1]        

        L = self.Lval
        
        Req = L / dt
        veq = -L * i[n - 1] / dt
        return {self.Reqsym:Req, self.Veqsym:veq}        


class SimulationResultsNode(object):

    def __init__(self, v):
        self.v = v


class SimulationResultsCpt(object):

    def __init__(self, v, i):
        self.v = v
        self.i = i

        
class SimulationResults(object):

    def __init__(self, tv, cct, r_model, node_list, branch_list):

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
            return SimulationResultsNode(self.node_voltages_get(name))
        
        if name in cct._elements:
            return SimulationResultsCpt(self.cpt_voltages_get(name),
                                        self.cpt_currents_get(name))

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

    def node_voltages_get(self, n):

        index = self.r_model._node_index(n)
        # NB, node_voltages is zero for index = -1
        return self.node_voltages[index]
        
    def cpt_voltages_get(self, cptname):

        cpt = self.cct.elements[cptname]        

        v1 = self.node_voltages_get(cpt.nodes[0])
        v2 = self.node_voltages_get(cpt.nodes[1])        
        return v1 - v2

    def cpt_voltage_get(self, cptname, n):

        cpt = self.cct.elements[cptname]
        
        v1 = self.node_voltages_get(cpt.nodes[0])[n]
        v2 = self.node_voltages_get(cpt.nodes[1])[n]     
        return v1 - v2

    def cpt_currents_get(self, cptname):

        try:
            index = self.r_model._branch_index(cptname)
            return self.branch_currents[index]
        except:
            cpt = self.cct._elements[cptname]
            if cpt.is_capacitor or cpt.is_inductor:
                # For a capacitor we can find the current through the
                # companion resistor or voltage source.
                return self.cpt_currents_get('V%seq' % cptname)

            Vd = self.cpt_voltages_get(cptname)            
            if cpt.is_resistor:
                return Vd / cpt.R.expr
                
            # Need to determine resistance of the cpt
            raise ValueError('FIXME')            

    def cpt_current_get(self, cptname, n):

        return self.cpt_currents_get(cptname)[n]
        
    @property
    def V(self, node):
        """Node voltage with respect to ground."""

        return self.get_Vd(node, '0')


    def get_Vd(self, Np, Nm=None):
        """Voltage drop between nodes"""

        self._solve()
        return (self._Vdict[Np] - self._Vdict[Nm]).canonical()    

    
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

        # Substitute values into the MNA A matrix and Z vector,
        # then perform numerical inversion of the A matrix.
        # Alternatively, a faster way, if the A matrix was small,
        # would be to symbolically invert the A matrix and
        # then substitute values.

        if n == 0:
            return

        dt = tv[n] - tv[n - 1]

        subsdict = {tsym: tv[n]}
        
        for cpt in self.reactive_cpts:

            # NB, node_voltages is zero for index = -1            
            v1 = results.node_voltages[cpt.v1_index]
            v2 = results.node_voltages[cpt.v2_index]            
            i = results.branch_currents[cpt.i_index]
            
            subsdict.update(cpt.subsdict(n, dt, v1, v2, i))

        A = self.A.subs(subsdict)
        Z = self.Z.subs(subsdict)

        if n == 1:
            symbols = A.free_symbols.union(Z.free_symbols)
            if symbols != set():
                raise ValueError('There are undefined symbols %s: use subs to replace with numerical values' % symbols)
        
        A1 = array(A).astype(float)
        Z1 = array(Z).astype(float)        
        
        A1inv = linalg.inv(A1)
        
        results1 = dot(A1inv, Z1).squeeze()

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

        if integrator == 'trapezoid':
            Ccls = SimulatedCapacitorTrapezoid
            Lcls = SimulatedInductorTrapezoid
        elif integrator == 'backward-euler':
            Ccls = SimulatedCapacitorBackwardEuler
            Lcls = SimulatedInductorBackwardEuler            
        else:
            raise ValueError('Unknown integrator ' + integrator)

        r_model = self.r_model

        # Construct MNA matrices.
        r_model._analyse()

        self.A = r_model._A
        self.Z = r_model._Z
        
        self.reactive_cpts = []        
        for key, elt in self.cct.elements.items():
            if not (elt.is_inductor or elt.is_capacitor):
                continue
            v1_index = r_model._node_index(elt.nodes[0])
            v2_index = r_model._node_index(elt.nodes[1])
            i_index = r_model._branch_index('V%seq' % elt.name)
            if elt.is_inductor:
                cls = Lcls
            else:
                cls = Ccls
            self.reactive_cpts.append(cls(elt, v1_index, v2_index, i_index))
        
        results = SimulationResults(tv, self.cct, r_model, r_model.node_list,
                                    r_model.unknown_branch_currents)
        
        for n, t1 in enumerate(tv):
            self._step(r_model, n, tv, results)

        return results

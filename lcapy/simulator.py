"""
This module implements an experimental time-stepping simulation.

Copyright 2020 Michael Hayes, UCECE
"""

from numpy import zeros, array, float, linalg, dot
from .sym import tsym, symbol_map


class SimulatedCapacitor(object):

    def __init__(self, C):

        self.nodes = C.nodes
        self.Cval = C.C.expr
        self.Reqname = 'R%seq' % C.name
        self.Veqname = 'V%seq' % C.name
        self.Reqsym = symbol_map(self.Reqname)
        self.Veqsym = symbol_map(self.Veqname)                

    
class SimulatedInductor(object):

    def __init__(self, L):

        self.nodes = L.nodes        
        self.Lval = L.L.expr
        self.Reqname = 'R%seq' % L.name
        self.Veqname = 'V%seq' % L.name
        self.Reqsym = symbol_map(self.Reqname)
        self.Veqsym = symbol_map(self.Veqname)        


class SimulatedCapacitorTrapezoid(SimulatedCapacitor):

    def subsdict(self, n, dt, results):
        """Create a dictionary of substitutions."""

        V1 = results.node_voltages_get(self.nodes[0])[n - 1]
        V2 = results.node_voltages_get(self.nodes[1])[n - 1]                

        V = V2 - V1
        C = self.Cval     

        Req = dt / (2 * C)
        Ieq = -2 * C * V / dt
        Veq = Req * Ieq

        return {self.Reqsym:Req, self.Veqsym:Veq}    


class SimulatedInductorTrapezoid(SimulatedInductor):

    def subsdict(self, n, dt, results):
        """Create a dictionary of substitutions."""        

        I = results.cpt_current_get(self.Veqname, n - 1)        
        L = self.Lval
        
        Req = 2 * L / dt
        Veq = -2 * L * I / dt
        return {self.Reqsym:Req, self.Veqsym:Veq}    


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
        
        self.node_voltages = zeros((self.num_nodes, N))
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
        if index < 0:
            return self.t * 0
        return self.node_voltages[index]
        
    def cpt_voltages_get(self, cptname):

        cpt = self.cct.elements[cptname]        

        V1 = self.node_voltages_get(cpt.nodes[0])
        V2 = self.node_voltages_get(cpt.nodes[1])        
        return V1 - V2

    def cpt_voltage_get(self, cptname, n):

        cpt = self.cct.elements[cptname]
        
        V1 = self.node_voltages_get(cpt.nodes[0])[n]
        V2 = self.node_voltages_get(cpt.nodes[1])[n]     
        return V1 - V2

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

    def __init__(self, cct, integrator='trapezoid'):

        if integrator == 'trapezoid':
            Ccls = SimulatedCapacitorTrapezoid
            Lcls = SimulatedInductorTrapezoid
        else:
            raise ValueError('Unknown integrator ' + integrator)
            
        self.cct = cct
        
        inductors = []
        capacitors = []

        for key, elt in cct.elements.items():
            if elt.is_inductor:
                inductors.append(elt)
            elif elt.is_capacitor:
                capacitors.append(elt)

        # Companion resistor model
        self.r_model = cct.r_model().subcircuits['time']

        self.capacitors = []
        for C1 in capacitors:
            self.capacitors.append(Ccls(C1))

        self.inductors = []
        for L1 in inductors:
            self.inductors.append(Lcls(L1))            

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
        
        for C1 in self.capacitors:
            subsdict.update(C1.subsdict(n, dt, results))

        for L1 in self.inductors:
            subsdict.update(L1.subsdict(n, dt, results))            

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
        results.node_voltages[:, n] = results1[0:num_nodes]
        results.branch_currents[:, n] = results1[num_nodes:]        

    def __call__(self, tv):

        N = len(tv)
        
        r_model = self.r_model

        # Construct MNA matrices.
        r_model._analyse()

        self.A = r_model._A
        self.Z = r_model._Z
        
        results = SimulationResults(tv, self.cct, r_model, r_model.node_list,
                                    r_model.unknown_branch_currents)
        
        for n, t1 in enumerate(tv):
            self._step(r_model, n, tv, results)

        return results
    

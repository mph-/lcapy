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
        #import pdb; pdb.set_trace()
        return {self.Reqsym:Req, self.Veqsym:Veq}    


class SimulatedInductorTrapezoid(SimulatedInductor):

    def subsdict(self, n, dt, results):
        """Create a dictionary of substitutions."""        

        I = results.cpt_current_get(self.Veqname, n - 1)        
        L = self.Lval
        
        Req = 2 * L / dt
        Veq = -2 * L * I / dt
        return {self.Reqsym:Req, self.Veqsym:Veq}    

    
class SimulationResults(object):

    def __init__(self, tv, cct, node_list, branch_list):

        self.t = tv
        self.cct = cct
        
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

        netlist = self.cct._netlist
        
        # If name is an integer, convert to a string.
        if isinstance(name, int):
            name = '%d' % name

        if name in netlist.nodes:
            return netlist.nodes[name]

        if name in netlist._elements:
            return netlist._elements[name]


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

        index = self.cct._node_index(n)
        if index < 0:
            return self.t * 0
        return self.node_voltages[index]
        
    def cpt_voltages_get(self, cptname):

        cpt = self.cct.elements[cptname]        

        V1 = self.node_voltages_get(cpt.nodes[0])
        V2 = self.node_voltages_get(cpt.nodes[1])        
        return V2 - V1

    def cpt_voltage_get(self, cptname, n):

        cpt = self.cct.elements[cptname]
        
        V1 = self.node_voltages_get(cpt.nodes[0])[n]
        V2 = self.node_voltages_get(cpt.nodes[1])[n]     
        return V2 - V1    

    def cpt_current_get(self, cptname, n):

        try:
            index = self.cct._branch_index(cptname)
            return self.branch_currents[index][n]
        except:
            Vd = self.voltage_get(cptname, n)
            # TODO: determine resistance of cpt but cannot do this for
            # a capacitor since it is time dependent.  So if desire a
            # capacitor current, will need to determine it as we go...
            # Alternatively, we need to record how the resistance
            # changes with time.  Or, we numerically approximate
            # i = C dv / dt.
            raise ValueError('FIXME')
        
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

        # Can symbolically invert the MNA matrix and perform
        # substitutions at each time step.  The symbolic matrix
        # inversion will be slow but the substitution is fast.  The
        # alternative is to substitute and then perform numerical
        # matrix inversion.  The latter is likely to be faster for
        # large problems.

        if n == 0:
            return

        dt = tv[n] - tv[n - 1]

        subsdict = {tsym: tv[n]}
        
        for C1 in self.capacitors:
            subsdict.update(C1.subsdict(n, dt, results))

        for L1 in self.inductors:
            subsdict.update(L1.subsdict(n, dt, results))            
            
        Z = self.Z.subs(subsdict)
        A = self.A.subs(subsdict)

        # TODO, need to evaluate functions such as Heaviside(t)

        # Could you lambdify but this is only advantageous if
        # evaluating a function for multiple values.
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
        
        results = SimulationResults(tv, r_model, r_model.node_list,
                                    r_model.unknown_branch_currents)
        
        for n, t1 in enumerate(tv):
            self._step(r_model, n, tv, results)

        return results
    

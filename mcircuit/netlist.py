"""
This module provides circuit analysis using modified nodal analysis
(MNA).

The circuit is described using netlists, similar to SPICE, but with
arbitrary node names (except for the ground node which is labelled 0).
The netlists can be loaded from a file or created at run-time.  For
example:

>>> from mcircuit import pprint, Circuit
>>> cct = Circuit('Voltage divider')
>>> cct.net_add('V_s fred 0') 
>>> cct.net_add('R_a fred 1') 
>>> cct.net_add('R_b 1 0') 
>>> cct.analyse()
>>> pprint(cct.V)
>>> pprint(cct.I)

cct.V is a directory of the nodal voltages.  cct.I is a directory of
the currents through independent voltage sources.  The directory keys
are the node names.  If the nodes are not integers, they need to
specified as strings:

>>> pprint(cct.V['fred'])
>>> pprint(cct.V[1])


Todo: 
1. Support dependent sources.



Copyright 2014 Michael Hayes, UCECE
"""

# SCApy  Symbolic Circuit Analysis in Python

from mcircuit import V, I, R, L, C, G, Is, Vs, pprint
import sympy as sym

# Implement modified nodal analysis (MNA)

class Element(object):


    def __init__(self, cpt, node1, node2, name):

        def node(name):
            
            # Convert to integer if possible
            try:
                node = int(name)
            except:
                node = name
            return node

        if not isinstance(cpt, (R, G, L, C, V, I)):
            raise ValueError('Adding component %s that is not R, G, L, C, V, I' % cpt)


        self.cpt = cpt
        self.name = name
        self.node1 = node(node1)
        self.node2 = node(node2)


    def __repr__(self):

        return 'Element(%s, %s, %s, %s)' % (self.cpt, self.node1, self.node2, self.name)


    def __str__(self):

        val = self.cpt.args[0]

        if val.is_symbol:
            return '%s %s %s' % (self.name, self.node1, self.node2)            

        return '%s %s %s %s' % (self.name, self.node1, self.node2, val.expr)


    @property
    def is_V(self):
        
        return isinstance(self.cpt, V)


    @property
    def is_I(self):
        
        return isinstance(self.cpt, I)


    @property
    def is_RLC(self):
        
        return isinstance(self.cpt, (R, G, L, C))


class NetElement(Element):

    def __init__(self, name, node1, node2, val=None):

        self.val = val
        if val == None:
            val = name

        kind = name[0]
        # An ammeter looks like a piece of wire so make a zero volt voltage source
        # so we can find the current through it.
        if kind == 'A':
            kind = 'V'
            val = 0
        
        # Allowable one-ports; this could be extended to Y, Z, etc.
        OPS = {'R' : R, 'G' : G, 'C' : C, 'L' : L, 'V' : V, 'I' : I}
        try:
            foo = OPS[kind]

        except KeyError:
            raise(ValueError, 'Unknown component kind for %s' % name)

        cpt = foo(val)

        super (NetElement, self).__init__(cpt, node1, node2, name)


    def __repr__(self):

        return 'NetElement(%s, %s, %s, %s)' % (self.name, self.node1, self.node2, self.val)


class Netlist(object):

    def __init__(self, filename=None):

        self.elements = []
        self.nodes = {}
        self.voltage_sources = []
        self.current_sources = []
        self.RLC = []

        self.V = {}
        self.I = {}
        self.cpt_counts = {'R' : 0, 'G' : 0, 'C' : 0, 'L' : 0, 'V' : 0, 'I' : 0}

        if filename != None:
            self.netfile_add(filename)


    def netfile_add(self, filename):    

        file = open(filename, 'r')
        
        lines = file.readlines()

        for line in lines:
            # Skip comments
            if line[0] in ('#', '%'):
                continue
            self.net_add(line.strip())


    def netlist(self):

        return '\n'.join([elt.__str__() for elt in self.elements])


    def _node_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = []
        self.nodes[node].append(elt)


    def _elt_add(self, elt):

        self.elements.append(elt)
        
        if elt.is_V: 
            self.voltage_sources.append(elt)
        if elt.is_I: 
            self.current_sources.append(elt)
        if elt.is_RLC: 
            self.RLC.append(elt)

        self._node_add(elt.node1, elt)
        self._node_add(elt.node2, elt)
        

    def net_add(self, line):

        parts = line.split(' ')
        
        elt = NetElement(*parts)

        self._elt_add(elt)


    def add(self, cpt, node1, node2, name=None):

        if name == None and hasattr(cpt.args[0], 'name'):
            name = cpt.args[0].name

        kind = type(cpt).__name__
        if not self.cpt_counts.has_key(kind):
            raise ValueError('Adding component %s that is not R, G, L, C, V, I' % cpt)

        if name == None:
            # Automatically generate unique name if one has not been specified
            self.cpt_counts[kind] = self.cpt_counts[kind] + 1
            name = '%s__%s' % (kind, self.cpt_counts[kind])

        elt = Element(cpt, node1, node2, name)

        self._elt_add(elt)
        
        
    def G_matrix_make(self):

        G = sym.zeros(self.num_nodes, self.num_nodes)

        for elt in self.RLC:
            n1 = elt.n1
            n2 = elt.n2
            Y = elt.cpt.Y

            if n1 >= 0 and n2 >= 0:
                G[n1, n2] -= Y
                G[n2, n1] -= Y
            if n1 >= 0:
                G[n1, n1] += Y
            if n2 >= 0:
                G[n2, n2] += Y

        return G


    def C_matrix_make(self):

        C = sym.zeros(len(self.voltage_sources), self.num_nodes)

        for m, elt in enumerate(self.voltage_sources):
            for n in range(self.num_nodes):
                if elt.n1 == n:
                    C[m, n] = 1
                elif elt.n2 == n:
                    C[m, n] = -1

        return C


    def D_matrix_make(self):

        D = sym.zeros(len(self.voltage_sources), len(self.voltage_sources))

        # Add dependent voltage sources here (VCVS and CCVS)
        return D


    def I_vector_make(self):

        I = sym.zeros(self.num_nodes, 1)

        for n in range(self.num_nodes):
            for m, elt in enumerate(self.current_sources):
                if elt.n1 == n:
                    I[n] = I[n] - elt.cpt.I
                elif elt.n2 == n:
                    I[n] = I[n] + elt.cpt.I
        return I


    def E_vector_make(self):

        E = sym.zeros(len(self.voltage_sources), 1)

        for m, elt in enumerate(self.voltage_sources):
            E[m] = elt.cpt.V
            
        return E


    def analyse(self):

        if not self.nodes.has_key(0):
            print('No ground node 0')

        self.nodemap = [0]
        self.revnodemap = {0 : 0}
        for n, node in enumerate(sorted(self.nodes.keys())):
            if node == 0:
                continue
            self.nodemap.append(node)
            self.revnodemap[node] = n

        self.num_nodes = len(self.nodemap) - 1

        # Assign mapped node numbers.  Note, ground is node -1.
        for elt in self.elements:
            elt.n1 = self.revnodemap[elt.node1] - 1
            elt.n2 = self.revnodemap[elt.node2] - 1

        # Compute matrices
        G = self.G_matrix_make()
        C = self.C_matrix_make()
        B = C.T
        D = self.D_matrix_make()

        # Augment the admittance matrix
        G = G.row_join(B).col_join(C.row_join(D))

        E = self.E_vector_make()
        I = self.I_vector_make()

        # Augment the known current vector with known voltage sources
        I = I.col_join(E)

        # Solve for the nodal voltages
        results = sym.simplify(G.inv() * I);        

        # Create dictionary of node voltages
        self.V = {}
        for n, node in enumerate(self.nodemap[1:]):        
            self.V[node] = Vs(results[n])

        # Create dictionary of currents through voltage sources
        self.I = {}
        for m, elt in enumerate(self.voltage_sources):
            self.I[elt.name] = Is(results[m + self.num_nodes])

        return self.V, self.I


class Circuit(Netlist):

    def __init__(self, circuitname):

        self.circuitname = circuitname
        super (Circuit, self).__init__()


def test():

    cct = Circuit('Test')

    cct.net_add('V_s fred 0') 
    cct.net_add('R_a fred bert') 
    cct.net_add('R_b bert 0') 
    cct.analyse()
    
    pprint(cct.V)
    
    pprint(cct.I)

    return cct



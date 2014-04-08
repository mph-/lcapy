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

from mcircuit import V, I, R, L, C, G, Vac, Iac, Is, Vs, pprint, cExpr
import sympy as sym

# Implement modified nodal analysis (MNA)

def _node(name):
            
    # Convert to integer if possible
    try:
        node = int(name)
    except:
        node = name
    return node


class CS(object):
    """Controlled source"""

    def __init__(self, node1, node2, arg):

        self.args = (node1, node2, arg)    
        # Controlling nodes
        self.cnodes = (_node(node1), _node(node2))
        self.arg = arg


class VCVS(CS):
    """Voltage controlled voltage source."""

    def __init__(self, node1, node2, A):
    
        A = cExpr(A)
        super (VCVS, self).__init__(node1, node2, A)
        # No independent component.
        self.V = 0


class Element(object):


    def __init__(self, cpt, node1, node2, name):

        if not isinstance(cpt, (R, G, L, C, V, I, Vac, Iac, VCVS)):
            raise ValueError('Adding component %s that is not R, G, L, C, V, I, VCVS' % cpt)


        self.cpt = cpt
        self.name = name
        self.nodes = (_node(node1), _node(node2))


    def __repr__(self):

        nodesstr = ', '.join(['%s' % node for node in self.nodes])
        return 'Element(%s, %s, %s)' % (self.cpt, nodesstr, self.name)


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes + self.args])


    @property
    def is_V(self):
        
        return isinstance(self.cpt, (V, Vac, VCVS))


    @property
    def is_I(self):
        
        return isinstance(self.cpt, (I, Iac))


    @property
    def is_RLC(self):
        
        return isinstance(self.cpt, (R, G, L, C))


class NetElement(Element):

    def __init__(self, name, node1, node2, *args):

        self.args = args

        if len(args) == 0:
            args = (name, )

        kind = name[0]
        # Handle Vac, Iac
        if len(name) > 3 and name[1:3] == 'ac':
            kind = name[0:3]

        # An ammeter looks like a piece of wire so make a zero volt voltage source
        # so we can find the current through it.
        if kind == 'A':
            kind = 'V'
            args = (0, )
        
        # Allowable one-ports; this could be extended to Y, Z, etc.
        OPS = {'R' : R, 'G' : G, 'C' : C, 'L' : L, 'V' : V, 'I' : I, 'Vac' : Vac, 'Iac' : Iac, 'E' : VCVS}
        try:
            foo = OPS[kind]

        except KeyError:
            raise(ValueError, 'Unknown component kind for %s' % name)

        cpt = foo(*args)

        super (NetElement, self).__init__(cpt, node1, node2, name)


    def __repr__(self):

        str = ', '.join(arg.__str__() for arg in [self.name] + list(self.nodes) + list(self.args))
        return 'NetElement(%s)' % str


class Netlist(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.nodes = {}
        self.voltage_sources = []
        self.current_sources = []
        self.RLC = []

        self.V = {}
        self.I = {}
        self.cpt_counts = {'R' : 0, 'G' : 0, 'C' : 0, 'L' : 0, 'V' : 0, 'I' : 0}

        if filename != None:
            self.netfile_add(filename)


    def __getitem__(self, key):

        return self.elements[key]


    def netfile_add(self, filename):    

        file = open(filename, 'r')
        
        lines = file.readlines()

        for line in lines:
            # Skip comments
            if line[0] in ('#', '%'):
                continue
            self.net_add(line.strip())


    def netlist(self):

        return '\n'.join([elt.__str__() for elt in self.elements.values()])


    def _node_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = []
        self.nodes[node].append(elt)


    def _elt_add(self, elt):

        if self.elements.has_key(elt.name):
            print('Overriding component %s' % elt.name)     
            # Need to search lists and update component.
           
        self.elements[elt.name] = elt

        self._node_add(elt.nodes[0], elt)
        self._node_add(elt.nodes[1], elt)
        

    def net_add(self, line):

        parts = line.split(' ')
        
        elt = NetElement(*parts)

        self._elt_add(elt)


    def add(self, cpt, node1, node2, name=None):
        """Add generic one-port component to netlist"""

        if name == None and hasattr(cpt.args[0], 'name'):
            name = cpt.args[0].name

        # Try to to collapse a composite component to a single component.
        cpt = cpt.simplify().cpt()

        kind = type(cpt).__name__
        if not self.cpt_counts.has_key(kind):
            raise ValueError('Adding unknown component %s' % cpt)

        if name == None:
            # Automatically generate unique name if one has not been specified
            self.cpt_counts[kind] = self.cpt_counts[kind] + 1
            name = '%s__%s' % (kind, self.cpt_counts[kind])

        elt = Element(cpt, node1, node2, name)

        self._elt_add(elt)
        
        
    def G_matrix(self):

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


    def B_matrix(self):

        C = sym.zeros(len(self.voltage_sources), self.num_nodes)

        for m, elt in enumerate(self.voltage_sources):
            for n in range(self.num_nodes):
                if elt.n1 == n:
                    C[m, n] = 1
                elif elt.n2 == n:
                    C[m, n] = -1

        return C.T


    def C_matrix(self):

        C = self.B_matrix().T

        for m, elt in enumerate(self.voltage_sources):
            if not isinstance(elt.cpt, VCVS):
                continue
            n1 = self.revnodemap[elt.cpt.cnodes[0]] - 1
            n2 = self.revnodemap[elt.cpt.cnodes[1]] - 1
            A = elt.cpt.arg
            
            for n in range(self.num_nodes):
                if n1 == n:
                    C[m, n] -= A
                elif n2 == n:
                    C[m, n] += A

        return C


    def D_matrix(self):

        D = sym.zeros(len(self.voltage_sources), len(self.voltage_sources))

        # Add dependent voltage sources here (VCVS and CCVS)
        return D


    def A_matrix(self):

        G = self.G_matrix()
        B = self.B_matrix()
        C = self.C_matrix()
        D = self.D_matrix()

        # Augment the admittance matrix to form A matrix
        A = G.row_join(B).col_join(C.row_join(D))

        return A


    def I_vector(self):

        I = sym.zeros(self.num_nodes, 1)

        for n in range(self.num_nodes):
            for m, elt in enumerate(self.current_sources):
                if elt.n1 == n:
                    I[n] = I[n] - elt.cpt.I
                elif elt.n2 == n:
                    I[n] = I[n] + elt.cpt.I
        return I


    def E_vector(self):

        E = sym.zeros(len(self.voltage_sources), 1)

        for m, elt in enumerate(self.voltage_sources):
            E[m] = elt.cpt.V
            
        return E


    def Z_vector(self):

        I = self.I_vector()
        E = self.E_vector()

        # Augment the known current vector with known voltage sources
        Z = I.col_join(E)
        return Z


    def analyse(self, mode='transient'):
        """mode either AC, DC, or transient"""

        if mode not in ('AC', 'DC', 'transient'):
            raise ValueError('Invalid analysis mode %s, must be AC, DC, transient' % mode)

        if not self.nodes.has_key(0):
            raise ValueError('Nothing connected to ground node 0')

        self.nodemap = [0]
        self.revnodemap = {0 : 0}
        for n, node in enumerate(sorted(self.nodes.keys())):
            if node == 0:
                continue
            self.nodemap.append(node)
            self.revnodemap[node] = n

        self.num_nodes = len(self.nodemap) - 1

        # Assign mapped node numbers.  Note, ground is node -1.
        for elt in self.elements.values():
            elt.n1 = self.revnodemap[elt.nodes[0]] - 1
            elt.n2 = self.revnodemap[elt.nodes[1]] - 1

        self.voltage_sources = []
        self.current_sources = []
        self.RLC = []
        for elt in self.elements.values():
            if elt.is_V: 
                self.voltage_sources.append(elt)
            elif elt.is_I: 
                self.current_sources.append(elt)
            elif elt.is_RLC: 
                if elt.cpt.V != 0.0:
                    print('Ignoring initial voltage on %s' % elt.name)
                    # Could use Norton model and split element into
                    # admittance in parallel with current source.
                self.RLC.append(elt)
            else:
                raise ValueError('Unhandled element %s' % elt.name)


        A = self.A_matrix()
        Z = self.Z_vector()

        # Solve for the nodal voltages
        results = sym.simplify(A.inv() * Z);        

        s, omega = sym.symbols('s omega')

        if mode in ('AC', 'DC'):
            results = results * s

        if mode == 'AC':
            # This a horrible interim hack...  Let's assume we are not
            # interested in the DC case; in this case DC voltage and current
            # sources zero at frequencies other than DC.
            # We should really only scale generic voltage or current sources
            # by s.
            results = results.subs(s, sym.I * omega)

        elif mode == 'DC':
            # This a better hack...  It will only work if all
            # voltage and current sources are DC.
            results = results.subs(s, sym.I * omega)
            results = results.subs(s, 0)

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



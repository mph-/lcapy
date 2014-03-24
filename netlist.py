from mcircuit import V, I, R, L, C, G
import sympy as sym

# Implement modified nodal analysis (MNA)

class Element(object):

    def __init__(self, name, node1, node2, val=None):

        self.kind = name[0]
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.val = val

        if val == None:
            val = self.name

        if self.kind == 'V':
            cpt = V(val)
        elif self.kind == 'I':
            cpt = I(val)
        elif self.kind == 'R':
            cpt = R(val)
        elif self.kind == 'G':
            cpt = G(val)
        elif self.kind == 'L':
            cpt = L(val)
        elif self.kind == 'C':
            cpt = C(val)
        else:
            raise(ValueError, 'Unknown component kind for %s' % name)

        self.cpt = cpt


    def __repr__(self):

        return 'Element(%s, %s, %s, %s)' % (self.name, self.node1, self.node2, self.val)


    @property
    def is_V(self):
        
        return isinstance(self.cpt, V)


    @property
    def is_I(self):
        
        return isinstance(self.cpt, I)


    @property
    def is_RLC(self):
        
        return isinstance(self.cpt, (R, G, L, C))


class Netlist(object):


    def __init__(self, filename):

        self.elements = []
        self.nodes = {}
        self.voltage_sources = []
        self.current_sources = []
        self.RLC = []

        file = open(filename, 'r')
        
        lines = file.readlines()

        for line in lines:
            
            parts = line.strip().split(' ')
            
            elt = Element(*parts)
            self.elements.append(elt)
            
            if elt.is_V: 
                self.voltage_sources.append(elt)
            if elt.is_I: 
                self.current_sources.append(elt)
            if elt.is_RLC: 
                self.RLC.append(elt)

        for elt in self.elements:
            self._node_add(elt.node1, elt)
            self._node_add(elt.node2, elt)

        
        if not self.nodes.has_key('0'):
            print('No ground node 0')

        self.nodemap = [0]
        self.revnodemap = {'0' : 0}
        for n, node in enumerate(sorted(self.nodes.keys())):
            if node == '0':
                continue
            self.nodemap.append(node)
            self.revnodemap[node] = n

        self.num_nodes = len(self.nodemap) - 1

        # Assign mapped node numbers
        for elt in self.elements:
            elt.n1 = self.revnodemap[elt.node1] - 1
            elt.n2 = self.revnodemap[elt.node2] - 1


    def _node_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = []
        self.nodes[node].append(elt)


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


    def I_vector_make(self):

        I = sym.zeros(self.num_nodes, 1)

        for n in range(self.num_nodes):
            for m, Is in enumerate(self.current_sources):
                if Is.n1 == n:
                    I[n] = I[n] - Is.cpt.I
                elif Is.n2 == n:
                    I[n] = I[n] + Is.cpt.I
        return I


    def V_vector_make(self):

        VV = sym.zeros(self.num_nodes, 1)

        for n, node in enumerate(self.nodemap[1:]):        
            VV[n] = V('V' + node)

        return VV


    def analyse(self):

        self.G = self.G_matrix_make()

        self.I = self.I_vector_make()
        self.V = self.V_vector_make()

        self.Vresult = sym.simplify(self.G.inv() * self.I);        

        Vresults = {}
        for n, node in enumerate(self.nodemap[1:]):        
            Vresults[node] = self.Vresult[n]
        return Vresults


def test():

    a = Netlist('net2.net')

    a.analyse()

    return a

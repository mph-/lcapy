"""
This module provides circuit analysis using modified nodal analysis
(MNA).

The circuit is described using netlists, similar to SPICE, but with
arbitrary node names (except for the ground node which is labelled 0).
The netlists can be loaded from a file or created at run-time.  For
example:

>>> from lcapy import pprint, Circuit
>>> cct = Circuit('Voltage divider')
>>> cct.net_add('V_s fred 0') 
>>> cct.net_add('R_a fred 1') 
>>> cct.net_add('R_b 1 0') 
>>> cct.analyse()
>>> pprint(cct.V)
>>> pprint(cct.I)

cct.V is a directory of the nodal voltages keyed by the node names.
If the nodes are not integers, they need to specified as strings.
cct.I is a directory of the currents through the components keyed by
the component names.  For example,

>>> pprint(cct.V['fred'])
>>> pprint(cct.V[1])
>>> pprint(cct.I['R1'])


Copyright 2014 Michael Hayes, UCECE
"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from warnings import warn
from lcapy.core import  pprint, cExpr
from lcapy.oneport import V, I, R, L, C, G, Y, Z, Vdc, Idc, Vac, Iac, Is, Vs, Ys, Zs
import sympy as sym


__all__ = ('Circuit', )


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


class TF(CS):
    """Ideal transformer.  T is turns ratio (secondary / primary)"""

    def __init__(self, node1, node2, T):
    
        T = cExpr(T)
        super (TF, self).__init__(node1, node2, T)
        # No independent component.
        self.V = 0


class Element(object):


    def __init__(self, cpt, node1, node2, name):

        if not isinstance(cpt, (R, G, L, C, V, I, Vac, Iac, VCVS, TF)):
            raise ValueError('Adding component %s that is not R, G, L, C, V, I, VCVS, TF' % cpt)


        self.cpt = cpt
        self.name = name
        self.nodes = (_node(node1), _node(node2))


    def __repr__(self):

        nodesstr = ', '.join(['%s' % node for node in self.nodes])
        return 'Element(%s, %s, %s)' % (self.cpt, nodesstr, self.name)


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes + self.cpt.args])


    @property
    def is_V(self):
        
        return isinstance(self.cpt, (V, Vdc, Vac, VCVS, TF))


    @property
    def is_I(self):
        
        return isinstance(self.cpt, (I, Idc, Iac))


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
        elif len(name) > 2 and name[0:2] == 'TF':
            kind = name[0:2]

        # An ammeter looks like a piece of wire so make a zero volt voltage source
        # so we can find the current through it.
        if kind == 'A':
            kind = 'V'
            args = (0, )
        
        # Allowable one-ports; this could be extended to Y, Z, etc.
        # Note V and I map to Vdc and Idc.
        OPS = {'R' : R, 'G' : G, 'C' : C, 'L' : L, 'V' : Vdc, 'I' : Idc, 'Vac' : Vac, 'Iac' : Iac, 'E' : VCVS, 'TF' : TF}
        try:
            foo = OPS[kind]

        except KeyError:
            raise(ValueError, 'Unknown component kind for %s' % name)

        cpt = foo(*args)

        super (NetElement, self).__init__(cpt, node1, node2, name)


    def __repr__(self):

        str = ', '.join(arg.__str__() for arg in [self.name] + list(self.nodes) + list(self.args))
        return 'NetElement(%s)' % str


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes + self.args])


class Netlist(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.nodes = {}
        self.num_nodes = 0
        self.voltage_sources = []
        self.current_sources = []
        self.RLC = []

        self._V = {}
        self._I = {}
        self.cpt_counts = {'R' : 0, 'G' : 0, 'C' : 0, 'L' : 0, 'V' : 0, 'I' : 0}

        if filename != None:
            self.netfile_add(filename)


    def __getitem__(self, name):
        """Return component by name"""

        return self.elements[name]


    def _nodeindex(self, node):
        """Return node index; ground is -1"""
        return self.revnodemap[node] - 1


    def netfile_add(self, filename):    
        """Add the nets from file with specified filename"""

        file = open(filename, 'r')
        
        lines = file.readlines()

        for line in lines:
            # Skip comments
            if line[0] in ('#', '%'):
                continue
            self.net_add(line.strip())


    def netlist(self):
        """Return the current netlist"""

        return '\n'.join([elt.__str__() for elt in self.elements.values()])


    def _node_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = []
        self.nodes[node].append(elt)


    def _invalidate(self):

        if hasattr(self, '_V'):
            delattr(self, '_V')
            delattr(self, '_I')


    def _elt_add(self, elt):

        if self.elements.has_key(elt.name):
            print('Overriding component %s' % elt.name)     
            # Need to search lists and update component.
           
        self._invalidate()

        self.elements[elt.name] = elt

        self._node_add(elt.nodes[0], elt)
        self._node_add(elt.nodes[1], elt)
        

    def net_add(self, line, *args):
        """The general form is: 'Name N1 N2 args'
        where N1 is the positive nose and N2 is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        parts = line.split(' ')
        
        elt = NetElement(*(tuple(parts) + args))

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
        

    def remove(self, name):
        """Remove specified element"""

        self._invalidate()

        if name not in self.elements:
            raise Error('Unknown component: ' + name)
        self.elements.pop(name)
        
        
    def G_matrix(self):

        G = sym.zeros(self.num_nodes, self.num_nodes)

        for elt in self.RLC:
            n1 = self._nodeindex(elt.nodes[0])
            n2 = self._nodeindex(elt.nodes[1])
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

        B = sym.zeros(self.num_nodes, len(self.voltage_sources))

        for m, elt in enumerate(self.voltage_sources):

            n1 = self._nodeindex(elt.nodes[0])
            n2 = self._nodeindex(elt.nodes[1])

            if n1 >= 0:
                B[n1, m] = 1
            if n2 >= 0:
                B[n2, m] = -1

            if isinstance(elt.cpt, TF):

                n3 = self._nodeindex(elt.cpt.cnodes[0])
                n4 = self._nodeindex(elt.cpt.cnodes[1])
                T = elt.cpt.arg
                
                if n3 >= 0:
                    B[n3, m] = -T
                if n4 >= 0:
                    B[n4, m] = T

        return B


    def C_matrix(self):

        C = sym.zeros(len(self.voltage_sources), self.num_nodes)

        for m, elt in enumerate(self.voltage_sources):

            n1 = self._nodeindex(elt.nodes[0])
            n2 = self._nodeindex(elt.nodes[1])

            if n1 >= 0:
                C[m, n1] = 1
            if n2 >= 0:
                C[m, n2] = -1

            if isinstance(elt.cpt, (TF, VCVS)):

                n3 = self._nodeindex(elt.cpt.cnodes[0])
                n4 = self._nodeindex(elt.cpt.cnodes[1])
                A = elt.cpt.arg
                
                if n3 >= 0:
                    C[m, n3] = -A
                if n4 >= 0:
                    C[m, n4] = A

        return C


    def D_matrix(self):

        D = sym.zeros(len(self.voltage_sources), len(self.voltage_sources))

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
                n1 = self._nodeindex(elt.nodes[0])
                n2 = self._nodeindex(elt.nodes[1])
                if n1 == n:
                    I[n] = I[n] - elt.cpt.I
                elif n2 == n:
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


    def _analyse(self):
        """Force reanalysis of network."""

        if not self.nodes.has_key(0):
            print('Nothing connected to ground node 0')
            self.nodes[0] = None

        self.nodemap = [0]
        self.revnodemap = {0 : 0}
        for n, node in enumerate(sorted(self.nodes.keys())):
            if node == 0:
                continue
            self.nodemap.append(node)
            self.revnodemap[node] = n

        self.num_nodes = len(self.nodemap) - 1

        self.voltage_sources = []
        self.current_sources = []
        self.RLC = []
        for elt in self.elements.values():
            if elt.is_V: 
                self.voltage_sources.append(elt)
            elif elt.is_I: 
                self.current_sources.append(elt)
            elif elt.is_RLC: 
                self.RLC.append(elt)
                if elt.cpt.V != 0.0:
                    # To handle initial condition, use Norton model
                    # and split element into admittance in parallel
                    # with current source.  The element will appear on
                    # both current_source and RLC lists.  We need to
                    # flip the current direction to follow convention
                    # that positive current flows from N1 to N2.
                    from lcapy import s

                    newelt = Element(I(elt.cpt.I), elt.nodes[1], elt.nodes[0],
                                     elt.name)
                    self.current_sources.append(newelt)
            else:
                raise ValueError('Unhandled element %s' % elt.name)


        A = self.A_matrix()
        Z = self.Z_vector()

        # Solve for the nodal voltages
        results = sym.simplify(A.inv() * Z);        

        # Create dictionary of node voltages
        self._V = {}
        self._V[0] = Vs(0)
        for n, node in enumerate(self.nodemap[1:]):        
            self._V[node] = Vs(results[n])

        # Create dictionary of currents through elements
        self._I = {}
        for m, elt in enumerate(self.voltage_sources):
            self._I[elt.name] = Is(results[m + self.num_nodes])

        for m, elt in enumerate(self.current_sources):
            self._I[elt.name] = elt.cpt.I

        # Don't worry about currents due to initial conditions; these
        # are overwritten below.
        for m, elt in enumerate(self.RLC):
            self._I[elt.name] = Is(sym.simplify((self._V[elt.nodes[0]] - self._V[elt.nodes[1]] - elt.cpt.V) / elt.cpt.Z))


    @property
    def V(self):    
        """Return dictionary of node voltages indexed by node name"""

        if not hasattr(self, '_V'):
            self._analyse()
        return self._V


    @property
    def I(self):    
        """Return dictionary of branch currents indexed by component name"""

        if not hasattr(self, '_I'):
            self._analyse()
        return self._I


    @property
    def Vd(self):    
        """Return dictionary of branch voltage differences indexed by component name"""
        if hasattr(self, '_Vd'):
            return self._Vd

        self._Vd = {}
        for elt in self.elements.values():
            self._Vd[elt.name] = Vs(sym.simplify(self.V[elt.nodes[0]] - self.V[elt.nodes[1]]))

        return self._Vd


    def Voc(self, n1, n2):
        """Determine open-circuit voltage between nodes."""

        return self.V[n1] - self.V[n2]

    
    def Isc(self, n1, n2):
        """Determine short-circuit current between nodes."""

        self.net_add('Vshort_', n1, n2, 0)

        Isc = self.I['Vshort_']
        self.remove('Vshort_')
        
        return Isc


    def thevenin(self, n1, n2):
        """Return Thevenin model between nodes n1 and n2"""

        Voc = self.Voc(n1, n2)
        Isc = self.Isc(n1, n2)

        return V(Voc) + Z(Zs(Voc / Isc))


    def norton(self, n1, n2):
        """Return Norton model between nodes n1 and n2"""

        Voc = self.Voc(n1, n2)
        Isc = self.Isc(n1, n2)

        return I(Isc) | Y(Ys(Isc / Voc))


    def Y(self, n1, n2):
        """Return admittance between nodes n1 and n2"""

        Voc = self.Voc(n1, n2)
        Isc = self.Isc(n1, n2)

        return Ys(Isc / Voc)


    def Z(self, n1, n2):
        """Return impedance between nodes n1 and n2"""

        Voc = self.Voc(n1, n2)
        Isc = self.Isc(n1, n2)

        return Zs(Voc / Isc)


class Circuit(Netlist):

    def __init__(self, circuitname=''):

        self.circuitname = circuitname
        super (Circuit, self).__init__()


def test():

    cct = Circuit('Test')

    cct.net_add('V_s fred 0') 
    cct.net_add('R_a fred bert') 
    cct.net_add('R_b bert 0') 
    
    pprint(cct.V)
    
    pprint(cct.I)

    return cct

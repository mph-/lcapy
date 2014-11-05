"""
This module provides circuit analysis using modified nodal analysis
(MNA).

The circuit is described using netlists, similar to SPICE, but with
arbitrary node names (except for the ground node which is labelled 0).
The netlists can be loaded from a file or created at run-time.  For
example:

>>> from lcapy import Circuit
>>> cct = Circuit('Voltage divider')
>>> cct.add('V_s fred 0') 
>>> cct.add('R_a fred 1') 
>>> cct.add('R_b 1 0') 
>>> cct.V.pprint()
>>> cct.I.pprint()

cct.V is a directory of the nodal voltages keyed by the node names.
If the nodes are not integers, they need to specified as strings.
cct.I is a directory of the currents through the components keyed by
the component names.  For example,

>>> cct.V['fred'].pprint()
>>> cct.V[1].pprint()
>>> cct.I['R1'].pprint()


Copyright 2014 Michael Hayes, UCECE
"""

# TODO: Add option to defer evaluation and thus keep things symbolic.
# This will help to simplify results that are not cancelled due to
# numerical quantisation.

from __future__ import division
from warnings import warn
from lcapy.core import  pprint, cExpr, Avs, Ais, Zs, Ys
from lcapy.oneport import V, I, v, i, R, L, C, G, Y, Z, Vdc, Idc, Vac, Iac, Is, Vs, Ys, Zs
from lcapy.twoport import AMatrix, TwoPortBModel
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
    def is_independentV(self):
        
        return isinstance(self.cpt, (V, Vdc, Vac))


    @property
    def is_independentI(self):
        
        return isinstance(self.cpt, (I, Idc, Iac))


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

        kind = name[0]
        if len(name) > 2 and name[0:2] == 'TF':
            kind = name[0:2]

        # Should check for Vdc1, etc.

        # Handle special cases for voltage and current sources.
        # Perhaps could generalise for impulse, step, ramp sources
        # although these can be specified symbolically, for example,
        # v1 1 0 t*Heaviside(t)
        if args != ():
            if kind == 'V' and args[0] == 'ac':
                kind = 'Vac'
            elif kind == 'V' and args[0] == 'dc':
                kind = 'Vdc'
            elif kind == 'I' and args[0] == 'ac':
                kind = 'Iac'
            elif kind == 'I' and args[0] == 'dc':
                kind = 'Idc'

            if kind in ('Vdc', 'Vac', 'Idc', 'Iac') and args[0] in ('ac', 'dc'):
                args = args[1:]

        # Allowable one-ports; this could be extended to Y, Z, etc.
        OPS = {'R' : R, 'G' : G, 'C' : C, 'L' : L, 'V' : V, 'Vdc' : Vdc, 'Vac' : Vac, 'v' : v, 'I' : I, 'Idc' : Idc, 'Iac' : Iac, 'i' : i, 'E' : VCVS, 'TF' : TF}
   
        try:
            foo = OPS[kind]

        except KeyError:
            raise(ValueError, 'Unknown component kind for %s' % name)

        if len(args) == 0:
            # Ensure symbol uppercase for s-domain value.
            if kind in ('Vdc', 'Vac', 'Idc', 'Iac'):
                name = name.capitalize()
            # Use component name for value
            args = (name, )

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
        self.voltage_sources = {}
        self.current_sources = {}
        self.RLC = []

        self._V = {0: Vs(0)}
        self._I = {}
        self.cpt_counts = {'R' : 0, 'G' : 0, 'C' : 0, 'L' : 0, 'V' : 0, 'I' : 0}

        if filename is not None:
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
            self.add(line.strip())


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
        

    def add(self, str, *args):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        parts = str.split(' ')
        elt = NetElement(*(tuple(parts) + args))
        self._elt_add(elt)


    def net_add(self, str, *args):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.

        Note, this method has been superseded by add.
        """

        parts = str.split(' ')
        elt = NetElement(*(tuple(parts) + args))
        self._elt_add(elt)


    def cpt_add(self, cpt, node1, node2, name=None):
        """Add generic one-port component to netlist"""

        if name is None and hasattr(cpt.args[0], 'name'):
            name = cpt.args[0].name

        # Try to to collapse a composite component to a single component.
        cpt = cpt.simplify().cpt()

        kind = type(cpt).__name__
        if not self.cpt_counts.has_key(kind):
            raise ValueError('Adding unknown component %s' % cpt)

        if name is None:
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
        
        if name in self.current_sources:
            self.current_sources.pop(name)

        if name in self.voltage_sources:
            self.voltage_sources.pop(name)

        
    def _G_matrix(self):

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


    def _B_matrix(self):

        B = sym.zeros(self.num_nodes, len(self.voltage_sources))

        for m, elt in enumerate(self.voltage_sources.values()):

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


    def _C_matrix(self):

        C = sym.zeros(len(self.voltage_sources), self.num_nodes)

        for m, elt in enumerate(self.voltage_sources.values()):

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


    def _D_matrix(self):

        D = sym.zeros(len(self.voltage_sources), len(self.voltage_sources))

        return D


    def _A_matrix(self):

        G = self._G_matrix()
        B = self._B_matrix()
        C = self._C_matrix()
        D = self._D_matrix()

        # Augment the admittance matrix to form A matrix
        A = G.row_join(B).col_join(C.row_join(D))

        return A


    def _I_vector(self):

        I = sym.zeros(self.num_nodes, 1)

        for n in range(self.num_nodes):
            for m, elt in enumerate(self.current_sources.values()):
                n1 = self._nodeindex(elt.nodes[0])
                n2 = self._nodeindex(elt.nodes[1])
                if n1 == n:
                    I[n] = I[n] - elt.cpt.I
                elif n2 == n:
                    I[n] = I[n] + elt.cpt.I
        return I


    def _E_vector(self):

        E = sym.zeros(len(self.voltage_sources), 1)

        for m, elt in enumerate(self.voltage_sources.values()):
            E[m] = elt.cpt.V
            
        return E


    def _Z_vector(self):

        I = self._I_vector()
        E = self._E_vector()

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

        self.voltage_sources = {}
        self.current_sources = {}
        self.RLC = []
        for key, elt in self.elements.iteritems():
            if elt.is_V: 
                self.voltage_sources[key] = elt
            elif elt.is_I: 
                self.current_sources[key] = elt
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
                    self.current_sources[key] = newelt
            else:
                raise ValueError('Unhandled element %s' % elt.name)


        A = self._A_matrix()
        Z = self._Z_vector()

        # Solve for the nodal voltages
        results = sym.simplify(A.inv() * Z);        

        # Create dictionary of node voltages
        self._V = {0: Vs(0)}
        for n, node in enumerate(self.nodemap[1:]):        
            self._V[node] = Vs(results[n])

        # Create dictionary of currents through elements
        self._I = {}
        for m, elt in enumerate(self.voltage_sources.values()):
            self._I[elt.name] = Is(results[m + self.num_nodes])

        for m, elt in enumerate(self.current_sources.values()):
            self._I[elt.name] = elt.cpt.I

        # Don't worry about currents due to initial conditions; these
        # are overwritten below.
        for m, elt in enumerate(self.RLC):
            self._I[elt.name] = Is(sym.simplify((self._V[elt.nodes[0]] - self._V[elt.nodes[1]] - elt.cpt.V) / elt.cpt.Z))


    @property
    def V(self):    
        """Return dictionary of s-domain node voltages indexed by node name"""

        if not hasattr(self, '_V'):
            self._analyse()
        return self._V


    @property
    def I(self):    
        """Return dictionary of s-domain branch currents indexed by component name"""

        if not hasattr(self, '_I'):
            self._analyse()
        return self._I


    @property
    def Vd(self):    
        """Return dictionary of s-domain branch voltage differences indexed by component name"""

        if hasattr(self, '_Vd'):
            return self._Vd

        self._Vd = {}
        for elt in self.elements.values():
            self._Vd[elt.name] = Vs(sym.simplify(self.V[elt.nodes[0]] - self.V[elt.nodes[1]]))

        return self._Vd


    @property
    def v(self):    
        """Return dictionary of t-domain node voltages indexed by node name"""

        if not hasattr(self, '_v'):
            self._v = {}
            for key in self.V:
                self._v[key] = self.V[key].inverse_laplace()

        return self._v


    @property
    def i(self):    
        """Return dictionary of t-domain branch currents indexed by component name"""

        if not hasattr(self, '_i'):
            self._i = {}
            for key in self.I:
                self._i[key] = self.I[key].inverse_laplace()

        return self._i


    @property
    def vd(self):    
        """Return dictionary of t-domain branch voltage differences indexed by component name"""

        if not hasattr(self, '_vd'):
            self._vd = {}
            for key in self.Vd:
                self._vd[key] = self.Vd[key].inverse_laplace()

        return self._vd


    def Voc(self, n1, n2):
        """Return open-circuit s-domain voltage between nodes n1 and n2."""

        return self.V[n1] - self.V[n2]


    def voc(self, n1, n2):
        """Return open-circuit t-domain voltage between nodes n1 and n2."""

        return self.Voc(n1, n2).inverse_laplace()

    
    def Isc(self, n1, n2):
        """Return short-circuit s-domain current between nodes n1 and n2."""

        self.add('Vshort_', n1, n2, 0)

        Isc = self.I['Vshort_']
        self.remove('Vshort_')
        
        return Isc


    def isc(self, n1, n2):
        """Return short-circuit t-domain current between nodes n1 and n2."""

        return self.Isc(n1, n2).inverse_laplace()


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
        """Return admittance between nodes n1 and n2
        with independent sources killed."""

        new = self.kill()

        Voc = new.Voc(n1, n2)
        Isc = new.Isc(n1, n2)

        return Ys(Isc / Voc)


    def Z(self, n1, n2):
        """Return impedance between nodes n1 and n2
        with independent sources killed."""

        new = self.kill()

        Voc = new.Voc(n1, n2)
        Isc = new.Isc(n1, n2)

        return Zs(Voc / Isc)


    def transfer(self, n1, n2, n3, n4):
        """Create voltage transfer function V2 / V1 where:
        V1 is V[n1] - V[n2]
        V2 is V[n3] - V[n4]
        
        Note, independent sources are killed."""

        new = self.kill()
        new.add('V1_', n1, n2)

        H = Avs(new.Voc(n3, n4) / new.Vd['V1_'])

        return H


    def Amatrix(self, n1, n2, n3, n4):
        """Create A matrix from network, where:
        I1 is the current flowing into n1 and out of n2
        I2 is the current flowing into n3 and out of n4
        V1 is V[n1] - V[n2]
        V2 is V[n3] - V[n4]
        """

        if self.Voc(n1, n2) != 0 or self.Voc(n3, n4) != 0:
            raise ValueError('Network contains independent sources')

        try:

            self.add('V1_', n1, n2)
            
            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = Avs(self.Vd['V1_'] / self.Voc(n3, n4))
            
            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = Zs(self.Vd['V1_'] / self.Isc(n3, n4))
            
            self.remove('V1_')
            
            self.add('I1_', n1, n2)
            
            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure I2 with port 2 open-circuit
            A21 = Ys(-self.I['I1_'] / self.Voc(n3, n4))
            
            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure I2 with port 2 short-circuit
            A22 = Ais(-self.I['I1_'] / self.Isc(n3, n4))
            
            self.remove('I1_')
            return AMatrix(A11, A12, A21, A22)

        except ValueError:
            raise ValueError('Cannot create A matrix')


    def kill(self):
        """Return a new circuit with the independent sources killed;
        i.e., make the voltage sources short-circuits and the current
        sources open-circuits."""

        new = Circuit()

        for key, elt in self.elements.iteritems():
            if elt.is_independentI: 
                continue
            if elt.is_independentV: 
                elt = NetElement(elt.name, elt.nodes[0], elt.nodes[1], 0)
            new._elt_add(elt)

        return new


    def twoport(self, n1, n2, n3, n4):
        """Create twoport model from network, where:
        I1 is the current flowing into n1 and out of n2
        I2 is the current flowing into n3 and out of n4
        V1 is V[n1] - V[n2]
        V2 is V[n3] - V[n4]
        """        

        V2b = self.Voc(n3, n4)
        I2b = self.Isc(n3, n4)

        A = self.kill().Amatrix(n1, n2, n3, n4)

        return TwoPortBModel(A.B, V2b, I2b)


class Circuit(Netlist):

    def __init__(self, circuitname=''):

        self.circuitname = circuitname
        super (Circuit, self).__init__()


def test():

    cct = Circuit('Test')

    cct.add('V_s fred 0') 
    cct.add('R_a fred bert') 
    cct.add('R_b bert 0') 
    
    pprint(cct.V)
    
    pprint(cct.I)

    return cct

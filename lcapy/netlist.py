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
from lcapy.core import  pprint, cExpr, Avs, Ais, Zs, Ys, s
from lcapy.oneport import V, I, v, i, Vdc, Idc, Vac, Iac, Vstep, Istep, Vacstep, Iacstep
from lcapy.oneport import R, L, C, G, Y, Z, Is, Vs, Ys, Zs
from lcapy.twoport import AMatrix, TwoPortBModel, Matrix, Vector
from schematic import Schematic
import sympy as sym
import re


__all__ = ('Circuit', )


cpt_types = ['C', # Capacitor
             'D', # Diode (not supported)
             'E', # VCVS
             'F', # CCCS (not supported yet, can be handled by G)
             'G', # VCCS (not supported yet)
             'H', # CCVS (not supported yet, can be handled by E)
             'I', # Current
             'K', # Mutual inductance
             'L', # Inductor
             'P', # Port (open-circuit)
             'Q', # Transistor (not supported)
             'R', # Resistor
             'TF', # Ideal transformer (even works at DC!)
             'TP', # Two-port (not supported yet)
             'V', # Voltage
             'W', # Wire (short-circuit)
             'Y', # Admittance
             'Z', # Impedance
         ]

# Note, SPICE netlists are usually case insensitive
# Perhaps prefix mechanical components with M?  But this will make
# crappy component labels.
mech_cpt_types = ['d', # Dashpot (damper, resistance)  perhaps b?
                  'f', # Force source
                  'k', # Spring
                  'm', # Mass
                  'u', # Velocity source
              ]

class Mdict(dict):

    def __getitem__(self, key):

        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key

        return super (Mdict, self).__getitem__(key)


class Ldict(dict):
    """Lazy dictionary for inverse Laplace"""

    def __init__(self, Vdict):

        self.Vdict = Vdict
        self.vdict = {}


    def __getitem__(self, key):
        
        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key

        if not self.vdict.has_key(key) and self.Vdict.has_key(key):
            self.vdict[key] = self.Vdict[key].inverse_laplace()

        return self.vdict[key]


    def keys(self):

        return self.Vdict.keys()


class CS(object):
    """Controlled source"""

    def __init__(self, *args):

        self.args = args   


class VCVS(CS):
    """Voltage controlled voltage source."""

    def __init__(self, A):
    
        A = cExpr(A)
        super (VCVS, self).__init__(A)
        # No independent component.
        self.V = 0


class TF(CS):
    """Ideal transformer.  T is turns ratio (secondary / primary)"""

    def __init__(self, T):
    
        T = cExpr(T)
        super (TF, self).__init__(T)
        # No independent component.
        self.V = 0


class K(object):
    """Mutual inductance"""

    def __init__(self, k):

        k = cExpr(k)
        self.k = k


class TP(CS):
    """Two-port Z-network"""

    def __init__(self, kind, Z11, Z12, Z21, Z22):
    
        super (TP, self).__init__(kind, Z11, Z12, Z21, Z22)
        # No independent component.
        self.V = 0


cpt_type_map = {'R' : R, 'C' : C, 'L' : L, 'Z' : Z, 'Y' : Y,
                'Vac' : Vac, 'Vdc' : Vdc, 
                'Iac' : Iac, 'Idc' : Idc, 
                'Vacstep' : Vacstep, 'Vstep' : Vstep,
                'Iacstep' : Iacstep, 'Istep' : Istep, 
                'Vimpulse' : V, 'Iimpulse' : I, 
                'Vs' : V, 'Is' : I, 
                'V' : v, 'I' : i, 'v' : v, 'i' : i,
                'P' : 'open', 'W' : 'short', 
                'E' : VCVS, 'TF' : TF, 'TP' : TP, 'K' : K, 
                'opamp' : VCVS}


# Regular expression alternate matches stop with first match so need
# to have longer names first.
cpt_types.sort(lambda x, y: cmp(len(y), len(x)))
cpt_type_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(cpt_types))


class Opts(dict):

    def _parse(self, string):

        for part in string.split(','):
            part = part.strip()
            if part == '':
                continue

            if part in ('up', 'down', 'left', 'right'):
                self['dir'] = part
                continue

            fields = part.split('=')
            key = fields[0].strip()
            arg = fields[1].strip() if len(fields) > 1 else ''
            self[key] = arg


    def __init__(self, arg):

        if isinstance(arg, str):
            self._parse(arg)
            return

        for key, val in arg.iteritems():
            self[key] = val


    def format(self):

        return ', '.join(['%s=%s' % (key, val) for key, val in self.iteritems()])


class Node(object):

    def __init__(self, name):

        self.name = name
        self.pos = None
        self.port = False
        parts = name.split('_')
        self.rootname = parts[0]  if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.list = []

    
    def append(self, elt):

        if elt.cpt_type in ('P', ):
            self.port = True

        self.list.append(elt)


class NetElement(object):

    def __init__(self, name, node1, node2, *args, **opts):

        match = cpt_type_pattern.match(name)

        if not match:
            raise ValueError('Unknown schematic component %s' % name)

        # Circuitikz does not like a . in a name
        if node1.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % node1)
        if node2.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % node2)

        cpt_type = match.groups()[0]
        id = match.groups()[1]

       
        self.opts = Opts(opts)
        self.name = name
        self.nodes = (node1, node2)
        self.args = args

        # Should check for Vdc1, etc.

        # Handle special cases for voltage and current sources.
        # Perhaps could generalise for impulse, step, ramp sources
        # although these can be specified symbolically, for example,
        # v1 1 0 t*Heaviside(t)


        if cpt_type == 'TP' and len(args) != 5:
            raise ValueError('TP component requires 5 args')

        cpt_type_orig = cpt_type
        if args != ():
            if cpt_type in ('V', 'I') and args[0] in ('ac', 'dc', 'step', 'acstep', 'impulse', 's'):
                cpt_type = cpt_type + args[0]
                args = args[1:]
            elif cpt_type == 'E' and args[0] == 'opamp':
                cpt_type = 'opamp'
                args = args[1:]

        if cpt_type in ('E', 'F', 'G', 'H', 'TF', 'TP', 'opamp'):
            if len(args) < 2:
                raise ValueError('Component type %s requires 4 nodes' % cpt_type)
            self.nodes += (args[0], args[1])
            args = args[2:]

        self.cpt_type = cpt_type

        if cpt_type in ('P', 'W'):
            self.cpt = None
            return

        try:
            foo = cpt_type_map[cpt_type]

        except KeyError:
            raise(ValueError, 'Unknown component kind for %s' % name)

        if len(args) == 0:
            # Ensure symbol uppercase for s-domain value.
            if cpt_type in ('Vdc', 'Vac', 'Idc', 'Iac'):
                name = name.capitalize()
            # Use component name for value
            args = (name, )

        cpt = foo(*args)
        self.cpt = cpt


    def __repr__(self):

        str = ', '.join(arg.__str__() for arg in [self.name] + list(self.nodes) + list(self.args))
        return 'NetElement(%s)' % str


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes[0:2] + self.args])


    @property
    def is_dummy(self):
        
        return self.cpt_type in ('P', 'W')


    @property
    def is_independentV(self):
        
        return isinstance(self.cpt, (V, Vdc, Vac, Vstep, Vacstep))


    @property
    def is_independentI(self):
        
        return isinstance(self.cpt, (I, Idc, Iac, Istep, Iacstep))


    @property
    def is_V(self):
        
        return isinstance(self.cpt, (V, Vdc, Vac, Vstep, Vacstep, VCVS, TF))


    @property
    def is_I(self):
        
        return isinstance(self.cpt, (I, Idc, Iac, Istep, Iacstep))


    @property
    def is_RC(self):
        
        return isinstance(self.cpt, (R, G, C))


    @property
    def is_L(self):
        
        return isinstance(self.cpt, L)


    @property
    def is_K(self):
        
        return isinstance(self.cpt, K)



class Netlist(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.nodes = {}
        # Shared nodes (with same voltage)
        self.snodes = {}

        self._V = Mdict({'0': Vs(0)})
        self._I = {}

        if filename is not None:
            self.netfile_add(filename)


    def __getitem__(self, name):
        """Return component by name"""

        return self.elements[name]


    def _node_index(self, node):
        """Return node index; ground is -1"""
        return self.node_list.index(self.node_map[node]) - 1


    def _branch_index(self, cpt_name):

        index = self.unknown_branch_currents.index(cpt_name)
        if index < 0:
            raise ValueError ('Unknown component name %s for branch current' % cpt.name)
        return index


    def netfile_add(self, filename):    
        """Add the nets from file with specified filename"""

        file = open(filename, 'r')
        
        lines = file.readlines()

        for line in lines:
            line = line.strip()
            if line == '':
                continue

            # Skip comments
            if line[0] in ('#', '%'):
                continue

            self.add(line)


    def netlist(self, full=False):
        """Return the current netlist"""

        from copy import copy

        lines = ''
        for key, elt in self.elements.iteritems():
            newelt = copy(elt)
            
            if not full:
                newelt.nodes = tuple([self.node_map[node] for node in elt.nodes])
                if elt.is_dummy:
                    continue

            line = newelt.__str__()
            if full:
                optstr = newelt.opts.format()
                if optstr != '':
                    line += ' ; ' + optstr

            lines += line + '\n'

        return lines


    def _node_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = Node(node)
        self.nodes[node].append(elt)

        vnode = self.nodes[node].rootname

        if vnode not in self.snodes:
            self.snodes[vnode] = []

        if node not in self.snodes[vnode]:
            self.snodes[vnode].append(node)


    def _invalidate(self):

        for attr in ('_V', '_I', '_node_map', '_node_list', '_sch', '_A', '_Z'):
            if hasattr(self, attr):
                delattr(self, attr)


    def _elt_add(self, elt):

        if self.elements.has_key(elt.name):
            print('Overriding component %s' % elt.name)     
            # Need to search lists and update component.
           
        self._invalidate()

        self.elements[elt.name] = elt

        # Ignore nodes for mutual inductance.
        if elt.cpt_type == 'K':
            return

        self._node_add(elt.nodes[0], elt)
        self._node_add(elt.nodes[1], elt)
        

    def net_parse(self, string, *args):

        fields = string.split(';')
        string = fields[1].strip() if len(fields) > 1 else ''
        if string != '':
            self.hints = True

        opts = Opts(string)

        parts = tuple(re.split(r'[\s]+', fields[0].strip()))
        elt = NetElement(*(parts + args), **opts)
        return elt


    def add(self, string, *args):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        elt = self.net_parse(string, *args)

        self._elt_add(elt)


    def net_add(self, string, *args):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.

        Note, this method has been superseded by add.
        """

        self.add(string, *args)


    def remove(self, name):
        """Remove specified element"""

        self._invalidate()

        if name not in self.elements:
            raise Error('Unknown component: ' + name)
        self.elements.pop(name)
        

    def _make_node(self):
        """Create a dummy node"""        

        if not hasattr(self, '_node_counter'):
            self._node_counter = 0
        self._node_counter += 1
        return '_%d' % self._node_counter


    def _make_open(self, node1, node2, opts):
        """Create a dummy open-circuit"""

        if not hasattr(self, '_open_counter'):
            self._open_counter = 0
        self._open_counter += 1

        net = 'P#%d %s %s ; %s' % (self._open_counter, node1, node2, opts.format())

        return self.net_parse(net)


    def _make_short(self, node1, node2, opts):
        """Create a dummy short-circuit"""

        if not hasattr(self, '_short_counter'):
            self._short_counter = 0
        self._short_counter += 1

        net = 'W#%d %s %s ; %s' % (self._short_counter, node1, node2, opts.format())

        return self.net_parse(net)


    def _make_Z(self, node1, node2, value, opts):
        """Create a dummy impedance"""

        if not hasattr(self, '_Z_counter'):
            self._Z_counter = 0
        self._Z_counter += 1

        net = 'Z#%d %s %s %s; %s' % (self._Z_counter, node1, node2, value, opts.format())

        return self.net_parse(net)


    def _make_V(self, node1, node2, value, opts):
        """Create a dummy s-domain voltage source"""

        if not hasattr(self, '_V_counter'):
            self._V_counter = 0
        self._V_counter += 1

        net = 'V#%d %s %s s %s; %s' % (self._V_counter, node1, node2, value, opts.format())

        return self.net_parse(net)


    def _make_I(self, node1, node2, value, opts):
        """Create a dummy s-domain current source"""

        if not hasattr(self, '_I_counter'):
            self._I_counter = 0
        self._I_counter += 1

        net = 'I#%d %s %s s %s; %s' % (self._I_counter, node1, node2, value, opts.format())

        return self.net_parse(net)


    def _RC_stamp(self, elt):
        """Add stamp for resistor or capacitor"""

        # L's can also be added with this stamp but if have coupling
        # it is easier to generate stamp that requires branch current
        # through the L.
        n1 = self._node_index(elt.nodes[0])
        n2 = self._node_index(elt.nodes[1])
        Y = elt.cpt.Y

        if n1 >= 0 and n2 >= 0:
            self._G[n1, n2] -= Y
            self._G[n2, n1] -= Y
        if n1 >= 0:
            self._G[n1, n1] += Y
        if n2 >= 0:
            self._G[n2, n2] += Y

        self._Is[n1] += elt.cpt.I


    def _L_stamp(self, elt):
        """Add stamp for inductor"""
        
        # This formulation adds the inductor current to the unknowns

        n1 = self._node_index(elt.nodes[0])
        n2 = self._node_index(elt.nodes[1])
        m = self._branch_index(elt.name)

        if n1 >= 0:
            self._B[n1, m] = 1
            self._C[m, n1] = 1
        if n2 >= 0:
            self._B[n2, m] = -1
            self._C[m, n2] = -1

        self._D[m, m] += -elt.cpt.Z

        self._Es[m] += elt.cpt.V


    def _K_stamp(self, elt):
        """Add stamp for mutual inductance"""

        # This requires the inductor stamp to include the inductor current.

        L1 = elt.nodes[0]
        L2 = elt.nodes[1]
        # TODO: Add sqrt to Expr
        ZM = elt.cpt.k * s * sym.simplify(sym.sqrt((self.elements[L1].cpt.Z * self.elements[L2].cpt.Z / s**2).expr))

        m1 = self._branch_index(L1)
        m2 = self._branch_index(L2)

        self._D[m1, m2] += -ZM
        self._D[m2, m1] += -ZM


    def _V_stamp(self, elt):
        """Add stamp for voltage source (independent and dependent)"""

        n1 = self._node_index(elt.nodes[0])
        n2 = self._node_index(elt.nodes[1])
        m = self._branch_index(elt.name)

        if n1 >= 0:
            self._B[n1, m] += 1
            self._C[m, n1] += 1
        if n2 >= 0:
            self._B[n2, m] -= 1
            self._C[m, n2] -= 1

        if isinstance(elt.cpt, TF):

            n3 = self._node_index(elt.nodes[2])
            n4 = self._node_index(elt.nodes[3])
            T = elt.cpt.args[0]
                
            if n3 >= 0:
                self._B[n3, m] -= T
                self._C[m, n3] -= T
            if n4 >= 0:
                self._B[n4, m] += T
                self._C[m, n4] += T

        elif isinstance(elt.cpt, VCVS):

            n3 = self._node_index(elt.nodes[2])
            n4 = self._node_index(elt.nodes[3])
            A = elt.cpt.args[0]
                
            if n3 >= 0:
                self._C[m, n3] -= A
            if n4 >= 0:
                self._C[m, n4] += A

        # Add ?
        self._Es[m] += elt.cpt.V


    def _I_stamp(self, elt):
        """Add stamp for current source (independent and dependent)"""

        n1 = self._node_index(elt.nodes[0])
        n2 = self._node_index(elt.nodes[1])
        if n1 >= 0:
            self._Is[n1] -= elt.cpt.I
        if n2 >= 0:
            self._Is[n2] += elt.cpt.I


    def _analyse(self):
        """Force reanalysis of network."""

        # TODO: think this out.  When a circuit is converted
        # to a s-domain model we get Z (and perhaps Y) components.
        # We also loose the ability to determine the voltage
        # across a capacitor or inductor since they get split
        # into a Thevenin model and renamed.
        if hasattr(self, '_s_model'):
            raise RuntimeError('Cannot analyse s-domain model')

        if not self.nodes.has_key('0'):
            print('Nothing connected to ground node 0')
            self.nodes['0'] = None


        # Determine which branch currents are needed.
        self.unknown_branch_currents = []

        for key, elt in self.elements.iteritems():
            if elt.is_V or elt.is_L:
                self.unknown_branch_currents.append(key)


        # Generate stamps.
        num_nodes = len(self.node_list) - 1
        num_branches = len(self.unknown_branch_currents)

        self._G = sym.zeros(num_nodes, num_nodes)
        self._B = sym.zeros(num_nodes, num_branches)
        self._C = sym.zeros(num_branches, num_nodes)
        self._D = sym.zeros(num_branches, num_branches)

        self._Is = sym.zeros(num_nodes, 1)
        self._Es = sym.zeros(num_branches, 1)

        for elt in self.elements.values():
            if elt.is_V:
                self._V_stamp(elt)
            elif elt.is_I: 
                self._I_stamp(elt)
            elif elt.is_RC: 
                self._RC_stamp(elt)
            elif elt.is_L: 
                self._L_stamp(elt)
            elif elt.is_K: 
                self._K_stamp(elt)
            elif elt.cpt_type not in ('P', 'W'):
                raise ValueError('Unhandled element %s' % elt.name)


        # Augment the admittance matrix to form A matrix
        self._A = self._G.row_join(self._B).col_join(self._C.row_join(self._D))
        # Augment the known current vector with known voltage vector
        # to form Z vector
        self._Z = self._Is.col_join(self._Es)


    def _solve(self):
        """Solve network."""

        if not hasattr(self, '_A'):
            self._analyse()

        # Solve for the nodal voltages
        try:
            Ainv = self._A.inv()
        except ValueError:
            raise ValueError('The A matrix is not invertible; probably some nodes need connecting with high value resistors')

        results = sym.simplify(Ainv * self._Z)

        # Create dictionary of node voltages
        self._V = Mdict({'0': Vs(0)})
        for n in self.nodes:
            index = self._node_index(n)
            if index >= 0:
                self._V[n] = Vs(results[index])
            else:
                self._V[n] = Vs(0)

        num_nodes = len(self.node_list) - 1

        # Create dictionary of branch currents through elements
        self._I = {}
        for m, key in enumerate(self.unknown_branch_currents):
            self._I[key] = Is(results[m + num_nodes])

        # Calculate the branch currents.  These should be evaluated as
        # required.  
        for key, elt in self.elements.iteritems():
            if elt.is_RC: 
                n1, n2 = self.node_map[elt.nodes[0]], self.node_map[elt.nodes[1]]
                self._I[elt.name] = Is(sym.simplify((self._V[n1] - self._V[n2] - elt.cpt.V) / elt.cpt.Z))



    @property
    def A(self):    
        """Return A matrix for MNA"""

        if not hasattr(self, '_A'):
            self._analyse()
        return Matrix(self._A)


    @property
    def Z(self):    
        """Return Z vector for MNA"""

        if not hasattr(self, '_Z'):
            self._analyse()
        return Vector(self._Z)


    @property
    def X(self):
        """Return X vector (of unknowns) for MNA"""

        if not hasattr(self, '_Z'):
            self._analyse()

        V = ['V_' + node for node in self.node_list[1:]]
        I = ['I_' + branch for branch in self.unknown_branch_currents]
        return Vector(V + I)


    @property
    def V(self):    
        """Return dictionary of s-domain node voltages indexed by node name"""

        if not hasattr(self, '_V'):
            self._solve()
        return self._V


    @property
    def I(self):    
        """Return dictionary of s-domain branch currents indexed by component name"""

        if not hasattr(self, '_I'):
            self._solve()
        return self._I


    @property
    def Vd(self):    
        """Return dictionary of s-domain branch voltage differences indexed by component name"""

        if hasattr(self, '_Vd'):
            return self._Vd

        self._Vd = {}
        for elt in self.elements.values():
            if elt.is_K:
                continue
            n1, n2 = self.node_map[elt.nodes[0]], self.node_map[elt.nodes[1]]
            self._Vd[elt.name] = Vs(sym.simplify(self.V[n1] - self.V[n2]))

        return self._Vd


    @property
    def v(self):    
        """Return dictionary of t-domain node voltages indexed by node name"""

        if not hasattr(self, '_v'):
            self._v = Ldict(self.V)

        return self._v


    @property
    def i(self):    
        """Return dictionary of t-domain branch currents indexed by component name"""

        if not hasattr(self, '_i'):
            self._i = Ldict(self.I)

        return self._i


    @property
    def vd(self):    
        """Return dictionary of t-domain branch voltage differences indexed by component name"""

        if not hasattr(self, '_vd'):
            self._vd = Ldict(self.Vd)

        return self._vd


    def Voc(self, Np, Nm):
        """Return open-circuit s-domain voltage between nodes Np and Nm."""

        return self.V[Np] - self.V[Nm]


    def voc(self, Np, Nm):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).inverse_laplace()

    
    def Isc(self, Np, Nm):
        """Return short-circuit s-domain current between nodes Np and Nm."""

        self.add('Vshort_ %d %d' %(Np, Nm), 0)

        Isc = self.I['Vshort_']
        self.remove('Vshort_')
        
        return Isc


    def isc(self, Np, Nm):
        """Return short-circuit t-domain current between nodes Np and Nm."""

        return self.Isc(Np, Nm).inverse_laplace()


    def thevenin(self, Np, Nm):
        """Return Thevenin model between nodes Np and Nm"""

        Voc = self.Voc(Np, Nm)

        # Connect 1 A s-domain current source between nodes and
        # measure voltage.
        self.add('Iin_ %d %d s 1' %(Nm, Np))
        Vf = self.Voc(Np, Nm)
        self.remove('Iin_')

        return V(Voc) + Z(Zs(Vf - Voc))


    def norton(self, Np, Nm):
        """Return Norton model between nodes Np and Nm"""

        Isc = self.Isc(Np, Nm)

        # Connect 1 V s-domain voltage source between nodes and
        # measure current.
        self.add('Vin_ %d %d s 1' %(Nm, Np))
        If = -self.I['Vin_']
        self.remove('Vin_')

        return I(Isc) | Y(Ys(If - Isc))


    def admittance(self, Np, Nm):
        """Return admittance between nodes Np and Nm
        with independent sources killed."""

        new = self.kill()

        Voc = new.Voc(Np, Nm)
        Isc = new.Isc(Np, Nm)

        return Ys(Isc / Voc)


    def impedance(self, Np, Nm):
        """Return impedance between nodes Np and Nm
        with independent sources killed."""

        new = self.kill()

        Voc = new.Voc(Np, Nm)
        Isc = new.Isc(Np, Nm)

        return Zs(Voc / Isc)


    def transfer(self, N1p, N1m, N2p, N2m):
        """Create voltage transfer function V2 / V1 where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        
        Note, independent sources are killed."""

        new = self.kill()
        new.add('V1_ %d %d impulse' % (N1p, N1m))

        H = Avs(new.Voc(N2p, N2m) / new.Vd['V1_'])

        return H


    def Amatrix(self, N1p, N1m, N2p, N2m):
        """Create A matrix from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """

        if self.Voc(N1p, N1m) != 0 or self.Voc(N2p, N2m) != 0:
            raise ValueError('Network contains independent sources')

        try:

            self.add('V1_ %d %d impulse' % (N1p, N1m))
            
            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = Avs(self.Vd['V1_'] / self.Voc(N2p, N2m))
            
            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure I2 with port 2 short-circuit
            A12 = Zs(self.Vd['V1_'] / self.Isc(N2p, N2m))
            
            self.remove('V1_')
            
            self.add('I1_ %d %d impulse' % (N1p, N1m))
            
            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure I2 with port 2 open-circuit
            A21 = Ys(-self.I['I1_'] / self.Voc(N2p, N2m))
            
            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure I2 with port 2 short-circuit
            A22 = Ais(-self.I['I1_'] / self.Isc(N2p, N2m))
            
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
                elt = self._make_open(elt.nodes[0], elt.nodes[1], elt.opts)
            elif elt.is_independentV: 
                elt = self._make_short(elt.nodes[0], elt.nodes[1], elt.opts)
            new._elt_add(elt)

        return new


    def twoport(self, N1p, N1m, N2p, N2m):
        """Create twoport model from network, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        """        

        V2b = self.Voc(N2p, N2m)
        I2b = self.Isc(N2p, N2m)

        A = self.kill().Amatrix(N1p, N1m, N2p, N2m)

        return TwoPortBModel(A.B, V2b, I2b)


    @property
    def node_map(self):
        """Determine mapping of nodes to unique nodes of the same potential"""

        if hasattr(self, '_node_map'):
            return self._node_map

        lnodes = self.lnodes

        node_map = {}
        for node in self.nodes:

            node_map[node] = node
            for nodes in lnodes:
                if node in nodes:
                    # Use first of the linked nodes unless '0' in list
                    if '0' in nodes:
                        node_map[node] = '0'
                    else:
                        # Sort nodes so 8 before 8_1 etc.
                        node_map[node] = sorted(nodes)[0]
                    break

        if not node_map.has_key('0'):
            print('Nothing connected to ground node 0')
            node_map['0'] = '0'

        self._node_map = node_map
        return node_map


    @property
    def node_list(self):
        """Determine list of unique nodes"""

        if hasattr(self, '_node_list'):
            return self._node_list

        # Extract unique nodes.
        node_list = list(set(self.node_map.values()))
        # Ensure node '0' is first in the list.
        node_list.insert(0, node_list.pop(node_list.index('0')))

        self._node_list = node_list
        return node_list


    @property
    def lnodes(self):
        """Determine linked nodes (both implicitly and explicitly
        connected)"""

        from copy import deepcopy

        # Start with implicitly linked nodes.
        lnodes = deepcopy(self.snodes)

        # Then augment with nodes connected by wires.
        for m, elt in enumerate(self.elements.values()):
            if elt.cpt_type not in ('W', ):
                continue

            n1, n2 = elt.nodes

            for key1, nodes in lnodes.iteritems():
                if n1 in nodes:
                    break;

            for key2, nodes in lnodes.iteritems():
                if n2 in nodes:
                    break;
                    
            if key1 != key2:
                lnodes[key1].extend(lnodes.pop(key2))

        # Remove nodes that are not linked.
        pnodes = []
        for key, nodes in lnodes.iteritems():
            if len(nodes) > 1:
                pnodes.append(nodes)

        return pnodes


    @property
    def sch(self):

        if hasattr(self, '_sch'):
            return self._sch

        netlist = self.netlist(full=True)

        sch = Schematic()

        for net in netlist.split('\n')[0:-1]:
            sch.add(net)

        self._sch = sch
        return sch


    def pre_initial_model(self):
        """Generate circuit model for determining the pre-initial conditions."""

        from copy import copy

        new_cct = self.__class__()

        for key, elt in self.elements.iteritems():
            
            # Assume initial C voltage and L current is zero.
                
            if elt.cpt_type in ('V', 'I', 'Vac', 'Iac'):
                print('Cannot determine pre-initial condition for %s, assuming 0' % elt.name)

            # v and i should be evaluated to determine the value at 0 - eps.
            if elt.cpt_type in ('v', 'i'):
                print('Cannot determine pre-initial condition for %s, assuming 0' % elt.name)

            if elt.cpt_type in ('C', 'Istep', 'Iacstep', 'I', 'i',
                                'Iac', 'Iimpulse'):
                elt = self._make_open(elt.nodes[0], elt.nodes[1], elt.opts)
            elif elt.cpt_type in ('L', 'Vstep', 'Vacstep', 'V', 'v',
                                  'Vac', 'Vimpulse'):
                elt = self._make_short(elt.nodes[0], elt.nodes[1], elt.opts)
            new_cct._elt_add(elt)

        return new_cct


    def s_model(self):

        from copy import copy

        cct = Circuit()
        cct._s_model = True

        for key, elt in self.elements.iteritems():
            
            new_elt = copy(elt)

            cpt_type = elt.cpt_type

            if cpt_type in ('C', 'L', 'R'):
                new_elt = self._make_Z(elt.nodes[0], elt.nodes[1], elt.cpt.Z, elt.opts)
            elif cpt_type in ('V', 'Vdc', 'Vac', 'Vimpulse', 'Vstep', 'Vacstep'):
                new_elt = self._make_V(elt.nodes[0], elt.nodes[1], elt.cpt.V, elt.opts)
            elif cpt_type in ('I', 'Idc', 'Iac', 'Iimpulse', 'Istep', 'Iacstep'):
                new_elt = self._make_I(elt.nodes[0], elt.nodes[1], elt.cpt.I, elt.opts)


            if cpt_type in ('C', 'L', 'R') and elt.cpt.V != 0:

                    dummy_node = self._make_node()


                    velt = self._make_V(dummy_node, elt.nodes[1], elt.cpt.V, elt.opts)
                    new_elt.nodes = (elt.nodes[0], dummy_node)

                    # Strip voltage label.  TODO: show voltage label across
                    # both components.
                    for opt in ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<'):
                        if new_elt.opts.has_key(opt):
                            new_elt.opts.pop(opt)

                    # Strip voltage and current labels.
                    for opt in ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<',
                                'i', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<', 
                                'i>_', 'i<_', 'i>^', 'i<^'):
                        if velt.opts.has_key(opt):
                            velt.opts.pop(opt)
  
                    cct._elt_add(velt)

            cct._elt_add(new_elt)
        
        return cct


    def draw(self, draw_labels=True, draw_nodes=True, label_nodes=True,
             s_model=False, filename=None, args=None, scale=1, stretch=1,
             tex=False):

        cct = self
        if s_model:
            cct = cct.s_model()
            
        return cct.sch.draw(draw_labels=draw_labels, draw_nodes=draw_nodes, 
                            label_nodes=label_nodes,
                            filename=filename, args=args, 
                            scale=scale, stretch=stretch, tex=tex)


class Circuit(Netlist):

    def __init__(self, filename=None):

        super (Circuit, self).__init__(filename)


def test():

    cct = Circuit('Test')

    cct.add('V_s fred 0') 
    cct.add('R_a fred bert') 
    cct.add('R_b bert 0') 
    
    pprint(cct.V)
    
    pprint(cct.I)

    return cct

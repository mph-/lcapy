"""
This module defines the components for modified nodal analysis.  The components
are defined at the bottom of this file.

Copyright 2015, 2016 Michael Hayes, UCECE

"""

from __future__ import print_function
from lcapy.core import cExpr, Vs, Is, s, sqrt
from copy import copy
import lcapy
import inspect
import sys

module = sys.modules[__name__]


class Cpt(object):

    def __init__(self, cct, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.cct = cct
        self.type = cpt_type
        self.id = cpt_id

        if cpt_id == '' and cct is not None:
            if cpt_type not in cct.anon:
                cct.anon[cpt_type] = 0
            cct.anon[cpt_type] += 1
            cpt_id = '#%d' % cct.anon[cpt_type]

        name = self.type + cpt_id

        self.string = string
        self.opts_string = opts_string
        self.nodes = nodes
        self.name = name
        self.args = args
        self.classname = self.__class__.__name__

        if self.type in ('W', 'O', 'P'):
            return

        if args is ():
            # Default value is the component name
            value = self.type
            if self.id != '':
                value += '_' + self.id

            args = (value, )
            self.args = args

        try:
            newclass = getattr(lcapy.oneport, self.classname)
        except:
            try:
                newclass = getattr(lcapy.twoport, self.classname)
            except:
                return
                
        self.cpt = newclass(*args)


    def __repr__(self):

        if hasattr(self, 'string'):
            return self.string
        
        return type(self)

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('stamp method not implemented for %s' % self)

    def kill_initial(self, newcct):
        """Copy cpt"""
        return newcct.add(self.string)

    def kill(self, newcct):
        raise ValueError('component not a source: %s' % self)        

    @property
    def I(self):
        """Current through element"""

        return self.cct.I[self.name]

    @property
    def i(self):
        """Time-domain current through element"""

        return self.cct.i[self.name]

    @property
    def V(self):
        """Voltage drop across element"""

        return self.cct.V[self.name]

    @property
    def v(self):
        """Time-domain voltage drop across element"""

        return self.cct.v[self.name]

    @property
    def Y(self):
        """Admittance"""
        
        return self.cpt.Y

    @property
    def Z(self):
        """Impedance"""
        
        return self.cpt.Z


    @property
    def node_indexes(self):

        return (self.cct._node_index(n) for n in self.nodes)

    @property
    def branch_index(self):

        return self.cct._branch_index(self.name)


class NonLinear(Cpt):

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('cannot analyse non-linear component %s' % self)


class TimeVarying(Cpt):

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('cannot analyse time-varying component %s' % self)


class O(Cpt):
    """Open circuit"""

    def stamp(self, cct, **kwargs):
        pass

class P(O):
    """Port"""
    pass


class RC(Cpt):

    def stamp(self, cct, **kwargs):

        # L's can also be added with this stamp but if have coupling
        # it is easier to generate stamp that requires branch current
        # through the L.
        n1, n2 = self.node_indexes

        Y = self.cpt.Y.expr

        if n1 >= 0 and n2 >= 0:
            cct._G[n1, n2] -= Y
            cct._G[n2, n1] -= Y
        if n1 >= 0:
            cct._G[n1, n1] += Y
        if n2 >= 0:
            cct._G[n2, n2] += Y

        if n1 >= 0:
            cct._Is[n1] += self.cpt.I.expr


class R(RC):
    pass


class C(RC):
    
    def kill_initial(self, newcct):
        # Kill implicit voltage sources due to initial conditions.
        return newcct.add('%s %s %s %s; %s' % (
            self.name, self.nodes[0], self.nodes[1], self.args[0], self.opts.format()))


class L(Cpt):
    
    def kill_initial(self, newcct):
        # Kill implicit voltage sources due to initial conditions.
        return newcct.add('%s %s %s %s; %s' % (
            self.name, self.nodes[0], self.nodes[1], self.args[0], self.opts.format()))

    def stamp(self, cct, **kwargs):

        # This formulation adds the inductor current to the unknowns

        n1, n2 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] = 1
            cct._C[m, n1] = 1
        if n2 >= 0:
            cct._B[n2, m] = -1
            cct._C[m, n2] = -1

        cct._D[m, m] += -self.cpt.Z.expr
        cct._Es[m] += self.cpt.V.expr


class E(Cpt):
    """VCVS"""

    def stamp(self, cct, **kwargs):
        n1, n2, n3, n4 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1
        
        A = cExpr(self.args[0]).expr
        
        if n3 >= 0:
            cct._C[m, n3] -= A
        if n4 >= 0:
            cct._C[m, n4] += A


class F(Cpt):
    """CCCS"""

    def stamp(self, cct, **kwargs):
        n1, n2 = self.node_indexes
        m = cct._branch_index(self.args[0])
        F = cExpr(self.args[1]).expr
            
        if n1 >= 0:
            cct._B[n1, m] -= F
        if n2 >= 0:
            cct._B[n2, m] += F


class G(Cpt):
    """VCCS"""

    def stamp(self, cct, **kwargs):
        n1, n2, n3, n4 = self.node_indexes
        G = cExpr(self.args[0]).expr

        if n1 >= 0 and n3 >= 0:
            cct._G[n1, n3] -= G
        if n1 >= 0 and n4 >= 0:
            cct._G[n1, n4] += G
        if n2 >= 0 and n3 >= 0:
            cct._G[n2, n3] += G
        if n2 >= 0 and n4 >= 0:
            cct._G[n2, n4] -= G


class H(Cpt):
    """CCVS"""

    def stamp(self, cct, **kwargs):
        n1, n2 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1
        
        mc = cct._branch_index(self.args[0])
        G = cExpr(self.args[1]).expr
        cct._D[m, mc] -= G


class I(Cpt):

    def kill(self, newcct):
        newopts = self.opts.copy()
        newopts.strip_voltage_labels()
        return newcct.add('O %s %s; %s' % (self.nodes[0], self.nodes[1], newopts.format()))


    def stamp(self, cct, **kwargs):

        n1, n2 = self.node_indexes

        I = self.cpt.I.expr
        if n1 >= 0:
            cct._Is[n1] += I
        if n2 >= 0:
            cct._Is[n2] -= I


class V(Cpt):

    def kill(self, newcct):
        newopts = self.opts.copy()
        newopts.strip_current_labels()
        return newcct.add('W %s %s; %s' % (self.nodes[0], self.nodes[1], newopts.format()))

    def stamp(self, cct, **kwargs):

        n1, n2 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1
        
        cct._Es[m] += self.cpt.V.expr


class K(Cpt):
    
    def __init__(self, cct, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super (K, self).__init__(cct, cpt_type, cpt_id, string, opts_string, nodes, *args)


    def stamp(self, cct, **kwargs):

        L1 = self.nodes[0]
        L2 = self.nodes[1]

        ZL1 = cct.elements[L1].cpt.Z
        ZL2 = cct.elements[L2].cpt.Z

        ZM = self.cpt.k * s * sqrt(ZL1 * ZL2 / s**2).simplify()

        m1 = cct._branch_index(L1)
        m2 = cct._branch_index(L2)

        cct._D[m1, m2] += -ZM.expr
        cct._D[m2, m1] += -ZM.expr


class TP(Cpt):
    """Two port"""
    pass


class TF(Cpt):
    """Transformer"""    

    def stamp(self, cct, **kwargs):

        n1, n2, n3, n4 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1
        
        T = self.cpt.args[0].expr

        if n3 >= 0:
            cct._B[n3, m] -= T
            cct._C[m, n3] -= T
        if n4 >= 0:
            cct._B[n4, m] += T
            cct._C[m, n4] += T


class W(Cpt):
    """Wire"""

    def stamp(self, cct, **kwargs):
        pass


class Y(RC):
    """Admittance"""
    pass


class Z(RC):
    """Impedance"""
    pass


classes = {}

def defcpt(name, base, docstring):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    classes[name] = newclass


def make(classname, parent, cpt_type, cpt_id,
         string, opts_string, nodes, *args):

    # Create instance of component object
    newclass = classes[classname]

    cpt = newclass(parent, cpt_type, cpt_id, string, opts_string, 
                   nodes, *args)
    # Add named attributes for the args?   Lname1, etc.
        
    return cpt


# Dynamically create classes.

defcpt('D', NonLinear, 'Diode')
defcpt('Dled', 'D', 'LED')
defcpt('Dphoto', 'D', 'Photo diode')
defcpt('Dschottky', 'D', 'Schottky diode')
defcpt('Dtunnel', 'D', 'Tunnel diode')
defcpt('Dzener', 'D', 'Zener diode')

defcpt('Eopamp', E, 'Opamp')
defcpt('Efdopamp', E, 'Fully differential opamp')

defcpt('sI', I, 's-domain current source')
defcpt('Isin', I, 'Sinusoidal current source')
defcpt('Idc', I, 'DC current source')
defcpt('Iac', I, 'AC current source')

defcpt('J', NonLinear, 'N JFET transistor')
defcpt('Jnjf', 'J', 'N JFET transistor')
defcpt('Jpjf', 'J', 'P JFET transistor')

defcpt('M', NonLinear, 'N MOSJFET transistor')
defcpt('Mnmos', 'M', 'N channel MOSJFET transistor')
defcpt('Mpmos', 'M', 'P channel MOSJFET transistor')

defcpt('Q', NonLinear, 'NPN transistor')
defcpt('Qpnp', 'Q', 'PNP transistor')
defcpt('Qnpn', 'Q', 'NPN transistor')

defcpt('SW', TimeVarying, 'Switch')
defcpt('SWno', 'SW', 'Normally open switch')
defcpt('SWnc', 'SW', 'Normally closed switch')
defcpt('SWpush', 'SW', 'Pushbutton switch')

defcpt('sV', V, 's-domain voltage source')
defcpt('Vsin', V, 'Sinusoidal voltage source')
defcpt('Vdc', V, 'DC voltage source')
defcpt('Vstep', V, 'Step voltage source')
defcpt('Vac', V, 'AC voltage source')

# Append classes defined in this module but not imported.
clsmembers = inspect.getmembers(module, lambda member: inspect.isclass(member) and member.__module__ == __name__)
for name, cls in clsmembers:
    classes[name] = cls


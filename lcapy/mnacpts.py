"""
This module defines the components for modified nodal analysis.  The components
are defined at the bottom of this file.

Copyright 2015, 2016 Michael Hayes, UCECE

"""

from __future__ import print_function
from lcapy.core import cExpr, s, sqrt, uppercase_name
from copy import copy
import lcapy
import inspect
import sys

module = sys.modules[__name__]


class Cpt(object):

    source = False

    def anon(self, cpt_type):

        cct = self.cct
        if cpt_type not in cct.anon:
            cct.anon[cpt_type] = 0
        cct.anon[cpt_type] += 1        
        return str(cct.anon[cpt_type])

    def __init__(self, cct, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.cct = cct
        self.type = cpt_type
        self.id = cpt_id

        if cpt_id == '' and cct is not None:
            cpt_id = '@' + self.anon(cpt_type)

        name = self.type + cpt_id

        self.net = string.split(';')[0]
        # This is the initial opts_string from which the opts attribute
        # is derived.
        self.opts_string = opts_string
        self.nodes = nodes
        self.name = name
        self.args = args
        self.classname = self.__class__.__name__
        self.opts = {}

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
        return self.__str__()

    def __str__(self):

        if self.opts == {}:
            return self.net
        return self.net + '; ' + str(self.opts)

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('stamp method not implemented for %s' % self)

    def kill_initial(self):
        """Kill implicit voltage sources due to initial conditions"""

        return str(self)

    def kill(self):
        """Kill component"""

        raise ValueError('component not a source: %s' % self)        

    def s_model(self, var):
        """Return s-domain model of component"""

        return str(self)

    def pre_initial_model(self):
        """Return pre-initial model of component"""

        return str(self)

    @property
    def causal(self):
        """Return True if causal component or if source produces
        a causal signal"""

        return self.cpt.causal

    @property
    def zeroic(self):
        """Return True if initial conditions are zero (or unspecified)"""

        return self.cpt.zeroic

    @property
    def hasic(self):
        """Return True if initial conditions are specified"""

        return self.cpt.hasic

    @property
    def I(self):
        """Current through component"""

        return self.cct.I[self.name]

    @property
    def i(self):
        """Time-domain current through component"""

        return self.cct.i[self.name]

    @property
    def V(self):
        """Voltage drop across component"""

        return self.cct.V[self.name]

    @property
    def v(self):
        """Time-domain voltage drop across component"""

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

    def dummy_node(self):

        return '_' + self.anon('node')


class NonLinear(Cpt):

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('cannot analyse non-linear component %s' % self)


class TimeVarying(Cpt):

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('cannot analyse time-varying component %s' % self)


class DummyCpt(Cpt):

    causal = True
    zeroic = True
    hasic = None


class O(DummyCpt):
    """Open circuit"""

    def stamp(self, cct, **kwargs):
        pass


class P(O):
    """Port"""
    pass


class RLC(Cpt):

    def s_model(self, var):

        if self.cpt.V == 0:        
            return 'Z%s %s %s {%s};%s' % (self.name, 
                                          self.nodes[0], self.nodes[1],
                                          self.cpt.Z(var), 
                                          self.opts)

        dummy_node = self.dummy_node()

        opts = self.opts.copy()

        # Strip voltage labels and save for open-circuit cpt
        # in parallel with Z and V.
        voltage_opts = opts.strip_voltage_labels()

        znet = 'Z%s %s %s {%s};%s' % (self.name, 
                                      self.nodes[0], dummy_node,
                                      self.cpt.Z(var), 
                                      opts)

        # Strip voltage and current labels from voltage source.
        opts.strip_all_labels()

        vnet = 'V%s %s %s s {%s}; %s' % (self.name, 
                                         dummy_node, self.nodes[1],
                                         self.cpt.V(var), opts)

        if voltage_opts == {}:
            return znet + '\n' + vnet

        # Create open circuit in parallel to the Z and V
        # that has the voltage labels.
        opts = self.opts.copy()
        opts.strip_current_labels()
        # Need to convert voltage labels to s-domain.
        # v(t) -> V(s)
        # v_C -> V_C
        # v_L(t) -> V_L(s)
        for opt, val in voltage_opts.items():
            opts[opt] = uppercase_name(val)
            
        onet = 'O%s %s %s; %s' % (self.name, 
                                  self.nodes[0], self.nodes[1], opts)
        return znet + '\n' + vnet + '\n' + onet


class RC(RLC):

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
    
    def kill_initial(self):
        """Kill implicit voltage sources due to initial conditions"""
        return '%s %s %s %s; %s' % (
            self.name, self.nodes[0], self.nodes[1], self.args[0], self.opts)

    def pre_initial_model(self):

        if self.cpt.v0 == 0.0:
            return 'O %s %s; %s' % (self.nodes[0], self.nodes[1], self.opts)
        return 'V%s %s %s %s; %s' % (self.name,
                                     self.nodes[0], self.nodes[1], 
                                     self.cpt.v0, self.opts)       


class L(RLC):
    
    def kill_initial(self):
        """Kill implicit voltage sources due to initial conditions"""
        return '%s %s %s %s; %s' % (
            self.name, self.nodes[0], self.nodes[1], self.args[0], self.opts)

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

    def pre_initial_model(self):

        if self.cpt.i0 == 0.0:
            return 'W %s %s; %s' % (self.nodes[0], self.nodes[1],
                                    self.opts)
        return 'I%s %s %s %s; %s' % (self.name,
                                     self.nodes[0], self.nodes[1], 
                                     self.cpt.i0, self.opts)       

class E(DummyCpt):
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


class F(DummyCpt):
    """CCCS"""

    def stamp(self, cct, **kwargs):
        n1, n2 = self.node_indexes
        m = cct._branch_index(self.args[0])
        F = cExpr(self.args[1]).expr
            
        if n1 >= 0:
            cct._B[n1, m] -= F
        if n2 >= 0:
            cct._B[n2, m] += F


class G(DummyCpt):
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


class H(DummyCpt):
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

    source = True

    def kill(self):
        newopts = self.opts.copy()
        newopts.strip_voltage_labels()
        return 'O %s %s; %s' % (self.nodes[0], self.nodes[1], newopts)

    def stamp(self, cct, **kwargs):

        n1, n2 = self.node_indexes

        I = self.cpt.I.expr
        if n1 >= 0:
            cct._Is[n1] += I
        if n2 >= 0:
            cct._Is[n2] -= I

    def s_model(self, var):

        return '%s %s %s s {%s}; %s' % (self.name, 
                                        self.nodes[0], self.nodes[1],
                                        self.cpt.I(var), self.opts)

    def pre_initial_model(self):

        # Assume IC zero.  FIXME
        return 'O %s %s; %s' % (self.nodes[0], self.nodes[1], self.opts)


class V(Cpt):

    source = True

    def kill(self):
        newopts = self.opts.copy()
        newopts.strip_current_labels()
        return 'W %s %s; %s' % (self.nodes[0], self.nodes[1], newopts)

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


    def s_model(self, var):

        return '%s %s %s s {%s}; %s' % (self.name, 
                                        self.nodes[0], self.nodes[1],
                                        self.cpt.V(var), self.opts)

    def pre_initial_model(self):

        # Assume IC zero.  FIXME
        return 'W %s %s; %s' % (self.nodes[0], self.nodes[1], self.opts)


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


class W(DummyCpt):
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

    # Switch context
    parent.context.switch()

    cpt = newclass(parent, cpt_type, cpt_id, string, opts_string, 
                   nodes, *args)
    # Add named attributes for the args?   Lname1, etc.

    # Restore context
    parent.context.restore()
        
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


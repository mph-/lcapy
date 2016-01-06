"""
This module defines the components for modified nodal analysis.  The components
are defined at the bottom of this file.

Copyright 2015, 2016 Michael Hayes, UCECE

"""

from __future__ import print_function
from lcapy.core import cExpr, Vs, Is, s, sqrt
from copy import copy

class Cpt(object):

    def __init__(self, sch, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.sch = sch
        self.type = cpt_type
        self.id = cpt_id

        if cpt_id is None:
            if cpt_type not in sch.anon:
                sch.anon[cpt_type] = 0
            sch.anon[cpt_type] += 1
            cpt_id = '#%d' % sch.anon[cpt_type]

        name = self.type + cpt_id

        self.string = string
        self.opts_string = opts_string
        self.nodes = nodes
        self.name = name
        self.args = args
        self.classname = self.__class__.__name__

    def __repr__(self):

        if hasattr(self, 'string'):
            return self.string
        
        return type(self)

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('stamp method not implemented for %s' % self)

    def kill(self):
        return copy(self)


class NonLinear(Cpt):

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('cannot analyse non-linear component %s' % self)


class TimeVarying(Cpt):

    def stamp(self, cct, **kwargs):
        raise NotImplementedError('cannot analyse time-varying component %s' % self)


class OpenCircuit(Cpt):

    def stamp(self, cct, **kwargs):
        pass


class RC(Cpt):

    def stamp(self, cct, **kwargs):

        # L's can also be added with this stamp but if have coupling
        # it is easier to generate stamp that requires branch current
        # through the L.
        n1 = cct._node_index(self.nodes[0])
        n2 = cct._node_index(self.nodes[1])

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


class Resistor(RC):
    pass


class Capacitor(RC):
    pass


class Inductor(Cpt):
    pass


class VCS(Cpt):
    pass


class CCS(Cpt):
    pass


class CurrentSource(Cpt):

    def kill(self):
        return Open()


class VoltageSource(Cpt):

    def kill(self):
        return Short()


class MutualInductance(Cpt):
    pass


class TwoPort(Cpt):
    pass


class Transformer(Cpt):
    pass


class Wire(Cpt):
    pass


class Admittance(Cpt):
    pass


class Impedance(Cpt):
    pass


classes = {}

def defcpt(name, base, docstring):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    classes[name] = newclass


# Dynamically create classes.

defcpt('C', Capacitor, 'Capacitor')

defcpt('D', NonLinear, 'Diode')
defcpt('Dled', 'D', 'LED')
defcpt('Dphoto', 'D', 'Photo diode')
defcpt('Dschottky', 'D', 'Schottky diode')
defcpt('Dtunnel', 'D', 'Tunnel diode')
defcpt('Dzener', 'D', 'Zener diode')

defcpt('E', VCS, 'VCVS')
defcpt('Eopamp', 'E', 'Opamp')
defcpt('Efdopamp', 'E', 'Fully differential opamp')
defcpt('F', VCS, 'VCCS')
defcpt('G', CCS, 'CCVS')
defcpt('H', CCS, 'CCCS')

defcpt('I', CurrentSource, 'Current source')
defcpt('Isin', 'I', 'Sinusoidal current source')
defcpt('Idc', 'I', 'DC current source')
defcpt('Iac', 'I', 'AC current source')

defcpt('J', NonLinear, 'N JFET transistor')
defcpt('Jnjf', 'J', 'N JFET transistor')
defcpt('Jpjf', 'J', 'P JFET transistor')

defcpt('K', MutualInductance, 'Mutual inductance')
defcpt('L', Inductor, 'Inductor')

defcpt('M', NonLinear, 'N MOSJFET transistor')
defcpt('Mnmos', 'M', 'N channel MOSJFET transistor')
defcpt('Mpmos', 'M', 'P channel MOSJFET transistor')

defcpt('O', OpenCircuit, 'Open circuit')
defcpt('P', OpenCircuit, 'Port')

defcpt('Q', NonLinear, 'NPN transistor')
defcpt('Qpnp', 'Q', 'PNP transistor')
defcpt('Qnpn', 'Q', 'NPN transistor')

defcpt('R', Resistor, 'Resistor')

defcpt('SW', TimeVarying, 'Switch')
defcpt('SWno', 'SW', 'Normally open switch')
defcpt('SWnc', 'SW', 'Normally closed switch')
defcpt('SWpush', 'SW', 'Pushbutton switch')

defcpt('TF', Transformer, 'Transformer')
defcpt('TP', TwoPort, 'Two port')

defcpt('V', VoltageSource, 'Voltage source')
defcpt('Vsin', 'V', 'Sinusoidal voltage source')
defcpt('Vdc', 'V', 'DC voltage source')
defcpt('Vstep', 'V', 'Step voltage source')
defcpt('Vac', 'V', 'AC voltage source')

defcpt('W', Wire, 'Wire')
defcpt('Y', Admittance, 'Admittance')
defcpt('Z', Impedance, 'Impedance')

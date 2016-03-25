"""This module provides circuit analysis using modified nodal analysis
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

Branch currents and branch voltage differences can be found using the
component name as an attribute, for example,

>>> cct.V_s.V.pprint()
>>> cct.R_a.I.pprint()

Nodal voltages (with respect to the ground node) can be found using
the node name or number as index, for example,

>>> cct['fred'].V.pprint()
>>> cct[1].V.pprint()

Copyright 2014, 2015, 2016 Michael Hayes, UCECE
"""

from lcapy.netlist import Netlist

__all__ = ('Circuit', )

class Circuit(Netlist):

    """Here's an example of using the Circuit class:

    cct = Circuit()
    cct.add('V1 1 0 V; down')
    cct.add('R1 1 2 R; right')
    cct.add('C1 2 0_2 C; down')
    cct.add('W 0 0_2; right')

    The directions after the semicolon are hints for drawing the
    schematic and are ignored for the circuit analysis.  The last net
    is a wire to make the schematic look nice; it is not needed for
    circuit analysis.  Indeed the capacitor could be connected
    directly to nodes 2 and 0.

    The nodes are usually numbers but can be any alphanumeric name
    including underscores.  By default, nodes with underscores are not
    drawn.

    The circuit components are also usually numbered but again they
    can be any alphanumeric name.  They can also have anonymous names,
    as for the wire in the example.  Internally they are enumerated
    sequentially for each component type: W#1, W#2, etc.

    The circuit can be displayed using:
    cct.draw()

    The schematic can be saved to a file using:
    cct.draw('schematic.pdf')

    The s-domain voltage across a component can be found using:
    cct.V1.V

    This is found using modified nodal analysis.  Once this is
    performed, the results are cached until the network is modified.

    The s-domain voltage through a component can be found using:
    cct.R1.I

    The s-domain nodal voltages with respect to the ground node (0)
    can be found using: 
    cct[2].V

    The time domain voltages and currents are displayed using
    lowercase attributes v and i.  For example,
    cct.C1.v

    Note that the answer assumes that all the dependent sources are
    zero for t < 0 and that all the inductors and capacitors have no
    initial currents and voltages.  Thus the Heaviside(t) factor
    should be ignored and replaced with the condition t >= 0.

    The impedance between nodes 2 and 0 can be found using:
    Z = cct.impedance(2, 0)

    The open-circuit voltage between nodes 2 and 0 can be found using:
    Z = cct.Voc(2, 0)

    The Thevenin equivalent circuit between nodes 2 and 0 can be found
    using:
    thevenin = cct.Thevenin(2, 0)

    The s-domain model can be drawn using:
    cct.s_model().draw()

    """

    def __init__(self, filename=None):

        super(Circuit, self).__init__(filename)

    def netfile_add(self, filename):
        """Add the nets from file with specified filename"""
        
        self._netfile_add(filename)

    def add(self, string):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        self._add(string)

    def Y(self, Np, Nm):
        """Return admittance between nodes Np and Nm with independent
        sources killed."""

        return self.admittance(Np, Nm)

    def Z(self, Np, Nm):
        """Return impedance between nodes Np and Nm with independent
        sources killed."""

        return self.impedance(Np, Nm)


def test():

    cct = Circuit('Test')

    cct.add('V_s fred 0')
    cct.add('R_a fred bert')
    cct.add('R_b bert 0')

    pprint(cct.V)

    pprint(cct.I)

    return cct

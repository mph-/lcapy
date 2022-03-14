"""This module provides circuit analysis using modified nodal analysis
(MNA).

The circuit is described using netlists, similar to SPICE, but with
arbitrary node names (except for the ground node which is labelled 0).
The netlists can be loaded from a file or created at run-time.  For
example:

>>> from lcapy import Circuit
>>> cct = Circuit('''
V_s fred 0
R_a fred 1
R_b 1 0''')

Branch currents and branch voltage differences can be found using the
component name as an attribute, for example,

>>> cct.V_s.V.pprint()
>>> cct.R_a.I.pprint()

Nodal voltages (with respect to the ground node) can be found using
the node name or number as index, for example,

>>> cct['fred'].V.pprint()
>>> cct[1].V.pprint()

Copyright 2014--2020 Michael Hayes, UCECE
"""

from .netlist import Netlist

__all__ = ('Circuit', )


class Circuit(Netlist):

    """The Circuit class is used for describing networks using
    netlists.  Despite the name, it does not require a closed path.

    Here's an example of using the Circuit class:

    cct = Circuit('''
    V1 1 0 V; down
    R1 1 2 R; right
    C1 2 0_2 C; down
    W 0 0_2; right''')

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

    The transform domain voltages across a component can be found using:
    cct.V1.V

    This is found using modified nodal analysis for each type of
    independent source in the circuit (AC, DC, transient, noise).
    Once this is performed, the results are cached until the network
    is modified.

    The transform domain currents through a component can be found using:
    cct.R1.I

    The transform domain nodal voltages with respect to the ground node (0)
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

    def __init__(self, filename=None, netlist=None, allow_anon=False, context=None):

        # Treat filename as netlist if it has a newline.
        if filename is not None and '\n' in filename:
            super(Circuit, self).__init__(allow_anon=allow_anon,
                                          context=context)
            self.add(filename)
        else:
            super(Circuit, self).__init__(filename, allow_anon=allow_anon,
                                          context=context)

        if netlist is not None:
            self.add(netlist)

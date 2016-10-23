.. _netlist:
.. _netlists:

========
Netlists
========

Circuits are described using a netlist of interconnected components (see :ref:`component-specification`).  Each line of a netlist describes a component using a Spice-like syntax.


Circuits
========

A circuit can be created by loading a netlist from a file or by
dynamically adding nets.  For example,

   >>> cct = Circuit('circuit.sch')

or

   >>> cct = Circuit()
   >>> cct.add('R1 1 2')


.. _component-specification:

Component specification
-----------------------

Each line in a netlist describes a single component, with the 
general form:

    component-name positive-node negative-node arg1 [arg2 etc.]

If no args are specified then the component value is assigned a
symbolic name specified by `component-name`. 

The component type is specified by the first letter(s) of the
`component-name`.  For example,

- DC voltage source of voltage V:

   Vname Np Nm dc V

- AC voltage source of voltage V, frequency f, and phase p:

   Vname Np Nm ac V f p

- Arbitrary s-domain voltage source:

   Vname Np Nm Vexpr

- DC current source of current I:

   Iname Np Nm dc I

- AC current source of current I, frequency f, and phase p:

   Iname Np Nm ac I p

- Arbitrary s-domain current source:

   Iname Np Nm Iexpr

- Resistor:

   Rname Np Nm R

- Conductor:

   Gname Np Nm G

- Inductor:

   Lname Np Nm L i0

- Capacitor:

   Cname Np Nm L v0

- Voltage-controlled voltage source (VCVS) of gain H with controlling nodes Nip and Nim:

   Ename Np Nm Nip Nim H

- Ideal transformer of turns ratio a:

   TFname Np Nm a

Np denotes the positive node; Np denotes the negative node.  Note,
positive current flows from `positive-node` to `negative-node`.  Node
names can be numeric or symbolic.  The ground node is designated `0`.


Voltage and current sources
---------------------------

The netlist description for a voltage source has the form:

Vname Np Nm value

Here value can be an arbitrary expression; the expression must be
enclosed in curly braces if it contains a delimiter such as a space,
comma, left bracket, or right bracket.  For example, for a ramp
voltage source

V1 1 0 {t * Heaviside(t)}

An s-domain value can be similarly described, for example

V1 1 0 {10 / s}

But what about the following example?

V1 1 2 10

Here the value is not an expression of t or s and so is ambiguous.
For Spice compatibility, Lcapy assumes a DC value of 10 (in the
s-domain this is equivalent to 10 / s).  An s-domain value of 10 can
be achieved using:

V1 1 2 {0 * s + 10}

or

V1 1 2 s 10

or alternatively,

V1 1 2 {10 * DiracDelta(t)}

To input an arbitrary time varying voltage use:

V1 1 2 {v(t)}

Here the value will be converted to V(s) for calculations but
displayed on a schematic as v(t).   

Here's an example of a cosine current current of amplitude 20 A and
frequency f

I1 1 0 {20 * cos(2 * pi * f * t)}

Here's an example of a negative exponential current of amplitude 20 A
and time constant 2 s

I1 1 0 {20 * exp(-t / 4)}


Netlist analysis examples
=========================


V-R-C circuit (1)
-----------------

This example plots the transient voltage across a capacitor in a series R-L circuit:

.. image:: examples/netlists/circuit-VRC1.png
   :width: 7cm

.. literalinclude:: examples/netlists/circuit-VRC1-vc.py

.. image:: examples/netlists/circuit-VRC1-vc.png
   :width: 15cm


V-R-C circuit (2)
-----------------

This example is the same as the previous example but it uses an
alternative method of plotting.

.. literalinclude:: examples/netlists/circuit-VRC2-vc.py

.. image:: examples/netlists/circuit-VRC2-vc.png
   :width: 15cm


V-R-L-C circuit (1)
-------------------

This example plots the transient voltage across a resistor in a series R-L-C circuit:

.. image:: examples/netlists/circuit-VRLC1.png
   :width: 7cm

.. literalinclude:: examples/netlists/circuit-VRLC1-vr.py

.. image:: examples/netlists/circuit-VRLC1-vr.png
   :width: 15cm



V-R-L-C circuit (2)
-------------------

This is the same as the previous example but with a different resistor value giving an underdamped response:

.. image:: examples/netlists/circuit-VRLC2-vr.png
   :width: 7cm

.. literalinclude:: examples/netlists/circuit-VRLC2-vr.py

.. image:: examples/netlists/circuit-VRLC2-vr.png
   :width: 15cm



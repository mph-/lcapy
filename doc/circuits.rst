================
Circuit Analysis
================


Introduction
============

Lcapy can only analyse linear time invariant (LTI) circuits, this
includes both passive and active circuits.

Time invariance means that the circuit parameters cannot change with
time; i.e., capacitors cannot change value with time.  It also means
that the circuit configuration cannot change with time, i.e., contain
switches (although switching problems can be analysed, see
:ref:`switching-analysis`).

Linearity means that superposition applies---if you double the voltage
of a source, the current (anywhere in the circuit) due to that source
will also double.  This restriction rules out components such as
diodes and transistors that have a non-linear relationship between
current and voltage (except in circumstances where the relationship
can be approximated as linear around some constant value---small
signal analysis).  Linearity also rules out capacitors where the
capacitance varies with voltage and inductors with hysteresis.

Superposition_ allows a circuit to be analysed by considering the
effect of each independent current and voltage source in isolation and
summing the results.


Linear circuit analysis
=======================

There is no universal analytical technique to determine the voltages
and currents in an LTI circuit.  Instead there are a number of methods
that all try to side step having to solve simultaneous
integro-differential equations.  These methods include DC analysis, AC
analysis, and Laplace analysis.


DC analysis
-----------

The simplest special case is for a DC independent source.  DC is an
idealised concept---it impossible to generate---but is a good enough
approximation for very slowly changing sources.  With a DC independent
source the dependent sources are also DC and thus no voltages or
currents change.  Thus capacitors can be replaced with open-circuits
and inductors can be replaced with short-circuits.  Currently, Lcapy
does not perform DC analysis although this can achieved to some extent
by using Lcapy to convert a circuit to a DC model.  Alternatively, a
Laplace analysis can be performed provided the initial conditions are
specified to avoid a transient reponse at :math:`t=0`.


AC analysis
-----------

With an AC source, phasor analysis can be performed.  Here the circuit
components are converted to an impedance dependent on the frequency of
the source. Currently, Lcapy does not perform AC analysis although
this can achieved to some extent by using Lcapy to convert a circuit
to an AC model.  Alternatively, a Laplace analysis can be performed
with :math:`s` replaced by :math:`\mathrm{j} \omega` assuming that the
initial conditions are specified to avoid a transient reponse at
:math:`t=0`.


Laplace analysis
----------------

The response due to a causal transient excitation from an independent
source can be analysed using Laplace analysis.  This is what Lcapy was
originally designed for.  A causal signal is zero for :math:`t<0`
analogous to a causal impulse response.

The response due to a non-causal excitation can be determined provided
the initial conditions (i.e., the voltages across capacitors and
currents through inductors) are specified.  More precisely, the
pre-intial conditions at :math:`t = 0_{-}` are required.  These differ
from the intial conditions at :math:`t = 0_{-}` whenever a Dirac delta
(or its derivative) excitation is considered.  Determining the initial
conditions is not straightforward for arbitrary excitations and at the
moment Lcapy expects you to do this!

The use of initial conditions also allows switching circuits to be
considered (see :ref:`switching-analysis`).

.. _switching-analysis:

Switching analysis
------------------

Whenever a circuit has a switch it is time variant.  The opening or
closing of switch changes the circuit and can produce transients.
While a switch violates the LTI requirements for linear circuit
analysis, the circuit prior to the switch changing can be analysed and
used to determine the initial conditions for the circuit after the
switched changed.  Lcapy requires that you do this!


Superposition
-------------

In principle, Lcapy could perform a combination of DC, AC, and Laplace
analysis to determine the overall result using superposition.
However, since the superposition needs to be performed in the time
domain, this requires a rethink of the s-domain attributes `V`, `I`,
`Y`, and `Z`.

Lcapy will happily kill of a specified independent source using the
`kill_except` method and thus s-domain superposition can be manually
performed.


.. _netlists:

Netlists
========

Circuits are described using a netlist of interconnected components (see :ref:`component-specification`).  Each line of a netlist describes a component using a Spice-like syntax.


.. _component-specification:

Component specification
-----------------------

The general form for a component is:

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

V1 1 2 {20 * DiracDelta(t)}

Unfortunately, the sympy Laplace transform does not consider the lower
limit of the integral to approach 0 in the limit from -oo and so the
DiracDelta gives half the value expected in circuit analysis.

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



Circuit examples
================


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



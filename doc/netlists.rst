.. _netlist:
.. _netlists:

========
Netlists
========

Circuits are described using a netlist of interconnected components (see :ref:`component-specification`).  Each line of a netlist describes a component using a Spice-like syntax.


Circuits
========

A circuit (or network) can be created by loading a netlist from a file or by
dynamically adding nets.  For example,

   >>> cct = Circuit('circuit.sch')

or

   >>> cct = Circuit()
   >>> cct.add('R1 1 2')
   >>> cct.add('L1 2 3')

or
   
   >>> cct = Circuit()
   >>> cct.add("""
   >>> R1 1 2
   >>> L1 2 3
   >>> """)

or
   
   >>> cct = Circuit("""
   >>> R1 1 2
   >>> L1 2 3
   >>> """)

This last version requires more than one net otherwise it is interpreted as a filename.   
   

.. _component-specification:

Component specification
-----------------------

Each line in a netlist describes a single component, with the 
general form:

    `component-name positive-node negative-node arg1 [arg2 etc.]`

If no args are specified then the component value is assigned a
symbolic name specified by `component-name`. 

Arguments containing delimiters (space, tab, comma, left bracket,
right bracket) can be escaped with brackets or double quotes.  For
example:

   `V1 1 0 {cos(5 * t)}`

The component type is specified by the first letter(s) of the
`component-name`.  For example,

- Arbitrary voltage source:

   `Vname Np Nm Vexpr`

   For example,

   `V1 1 0`  This is equivalent to `V1 1 0 {v1(t)}`

   `V1 1 0 10`  This is a DC source of 10 V
  
   `V1 1 0 {2 * cos(5 * t)}` This is an AC source
  
   `V1 1 0 {2 * cos(5 * t) * u(t)}`  This is a transient source

   `V1 1 0 {10 / s}` This is a transient source defined in the s-domain

   `V1 1 0 {s * 0 + 10}`  This is a transient source defined in the s-domain, equivalent to `V1 1 0 s 10`

- DC voltage source of voltage V:

   `Vname Np Nm dc V`

- AC voltage source of complex voltage amplitude V and phase p (radians) with angular frequency :math:`\omega`:  (omega)

   `Vname Np Nm ac V p`

- AC voltage source of complex voltage amplitude V and phase p (radians) with angular frequency w:

   `Vname Np Nm ac V p w`  

- Step voltage source of amplitude V

   `Vname Np Nm step V`

- s-domain voltage source of complex amplitude V

   `Vname Np Nm s V`  
   
- Arbitrary current source:

   `Iname Np Nm Iexpr`

   `I1 1 0`  This is equivalent to `I1 1 0 {v1(t)}`

   `I1 1 0 10`  This is a DC source of 10 I
  
   `I1 1 0 {2 * cos(5 * t)}` This is an AC source
  
   `I1 1 0 {2 * cos(5 * t) * u(t)}`  This is a transient source

   `I1 1 0 {10 / s}` This is a transient source defined in the s-domain

   `I1 1 0 {s * 0 + 10}`  This is a transient source defined in the s-domain, equivalent to `I1 1 0 s 10`

- DC current source of current I:

   `Iname Np Nm dc I`

- AC current source of complex current amplitude I and phase p (radians) with angular frequency :math:`\omega`:  (omega)

   `Iname Np Nm ac I p`

- AC current source of complex current amplitude I and phase p (radians) with angular frequency w:

   `Iname Np Nm ac I p w`

- Step current source of amplitude I

   `Iname Np Nm step I`

- s-domain current source of complex current I

   `Iname Np Nm s I`

- Resistor:

   `Rname Np Nm R`

- Conductor:

   `Gname Np Nm G`

- Inductor:

   `Lname Np Nm L`
  
   `Lname Np Nm L i0`  Here `i0` is the initial current through the inductor.  If this is specified then the circuit is solved as an initial value problem.

- Capacitor:

   `Cname Np Nm L`
 
   `Cname Np Nm L v0`   Here `i0` is the voltage across the capacitor.  If this is specified then the circuit is solved as an initial value problem.

- Voltage-controlled voltage source (VCVS) of gain H with controlling nodes Nip and Nim:

   `Ename Np Nm Nip Nim H`

- Ideal transformer of turns ratio a:

   `TFname Np Nm Nip Nim a`

- Ideal gyrator of gyration resistance R:

   `GYname Np Nm Nip Nim R`  

Np denotes the positive node; Np denotes the negative node.  For
two-port devices, Nip denotes the positive input node and Nim denotes
the negative input node.  Note, positive current flows from
`positive-node` to `negative-node`.  Node names can be numeric or
symbolic.  The ground node is designated `0`.

If the value is not explicity specified, the component name is used.
For example,

   `C1 1 0` is equivalent to `C1 1 0 C1`


Circuit attributes
------------------

A circuit is comprised of a collection of nodes and a collection of
circuit elements (components).  For example,

   >>> cct = Circuit()
   >>> cct.add('V1 1 0 {u(t)}')
   >>> cct.add('R1 1 2')
   >>> cct.add('L1 2 0')
   >>>
   >>> cct
   V1 1 0 {u(t)}
   R1 1 2
   L1 2 0

A node object is obtained using indexing notation, for example:

   >>> cct[2]

A circuit element object is obtained using its name, for example:

   >>> cct.R1   


Node attributes
---------------

Nodes have two attributes: `v`, `V`, `dpY`, and `dpZ`.

`v` is the time-domain voltage (with respect to the ground node 0).
`V` is a superposition of the node voltage in the different transform
domains.

For example,
   
   >>> cct[2].v
    -R₁⋅t              
    ──────             
      L₁               
   e      ⋅Heaviside(t)


`dpY` and `dpZ` return the driving-point admittance and impedance for the node with respect to ground, for example,

   >>> cct[2].dpZ
    R₁⋅s 
   ──────
       R₁
   s + ──
       L₁


   
Component attributes
--------------------

Circuit elements (components) have attributes: `v`, `V`, `i`, `I`,
`Y`, `Z`, `dpY`, and `dpZ`.

`v` is the time-domain voltage difference across the component, for example:

   >>> cct.R1.v   
    -R₁⋅t              
    ──────             
      L₁               
   e      ⋅Heaviside(t)   
   
`i` is the time-domain current through the component, for example:

   >>> cct.R1.i   
   ⎛      -R₁⋅t ⎞             
   ⎜      ──────⎟             
   ⎜        L₁  ⎟             
   ⎜1    e      ⎟             
   ⎜── - ───────⎟⋅Heaviside(t)
   ⎝R₁      R₁  ⎠             

The `V` and `I` attributes provide the voltage and current as a
superposition in the transform domains, for example,

   >>> cct.V1.V
   ⎧   1⎫
   ⎨s: ─⎬
   ⎩   s⎭

The `Y` and `Z` attributes provide the admittance and impedance of the
component, for example,

   >>> cct.L1.Z
   L₁⋅s

   >>> cct.R1.Z
   R₁

The driving point admittance and impedance can be found using `dpY` and `dpZ`, for example,

   >>> cct.L1.dpZ
    R₁⋅s 
   ──────
       R₁
   s + ──
       L₁

Note, this is the total impedance across `L1`.
       

   
Circuit evaluation
------------------

The circuit node voltages are determined using Modified Nodal Analysis
(MNA).  This is performed lazily as required with the results cached.

When a circuit has multiple independent sources, the circuit is
decomposed into a number of sub-circuits; one for each source type.
Again, this is performed lazily as required.  Each sub-circuit is
evaluated independently and the results are summed using the principle
of superposition.  For example, consider the circuit

   >>> cct = Circuit()
   >>> cct.add('V1 1 0 {1 + u(t)}')
   >>> cct.add('R1 1 2')
   >>> cct.add('L1 2 0')

In this example, V1 can be considered the superposition of a DC source
and a transient source.  The approach Lcapy uses to solve the circuit
can be found using the `describe` method:

   >>> cct.describe()
   This is solved using superposition.
   DC analysis is used for source V1.
   Laplace analysis is used for source V1.

For the curious, the sub-circuits can be found with the `sub` attribute:

   >>> cct.sub
   {'dc': V1 1 0 dc {1}
          R1 1 2
          L1 2 0 L_1,
   's': V1 1 0 {Heaviside(t)}
        R1 1 2
        L1 2 0 L_1
   }

Here the first sub-circuit is solved using DC analysis and the second
sub-circuit is solved using Laplace analysis in the s-domain.

The properties of each sub-circuit can be found with the `analysis` attribute:

   >>> cct.sub['dc'].analyse()
   {'ac': False,
   'causal': False,
   'control_sources': [],
   'dc': True,
   'dependent_sources': [],
   'has_s': False,
   'hasic': False,
   'independent_sources': ['V1'],
   'ivp': False,
   'time_domain': False,
   'zeroic': True}


Netlist analysis examples
=========================


V-R-C circuit (1)
-----------------

This example plots the transient voltage across a capacitor in a series R-L circuit:

.. image:: examples/netlists/circuit-VRC1.png
   :width: 10cm

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
   :width: 10cm

.. literalinclude:: examples/netlists/circuit-VRLC1-vr.py

.. image:: examples/netlists/circuit-VRLC1-vr.png
   :width: 15cm



V-R-L-C circuit (2)
-------------------

This is the same as the previous example but with a different resistor value giving an underdamped response:

.. image:: examples/netlists/circuit-VRLC2.png
   :width: 10cm
           
.. literalinclude:: examples/netlists/circuit-VRLC2-vr.py

.. image:: examples/netlists/circuit-VRLC2-vr.png
   :width: 15cm





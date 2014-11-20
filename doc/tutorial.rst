========
Tutorial
========


Introduction
============

lcapy is a Python package for linear circuit analysis.  It will only
solve linear, time invariant networks.  In other words, networks
comprised of basic circuit elements (R, L, C, etc.) that do not vary
with time.

It does not support non-linear devices such as diodes or transistors
although it does support simple opamps without saturation.

Lcapy uses Sympy (symbolic Python) for its values and expressions
and thus the circuit analysis can be performed symbolically.  See http://docs.sympy.org/latest/tutorial/index.html for the SymPy tutorial.

Internally, the circuit components are stored using their s-domain
equivalents, such as impedances and admittances.  This is convenient
for frequency response analysis but requires an inverse Laplace
transform for transient response analysis.


Preliminaries
=============

- Before you can use Lcapy you need to install the Lcapy package (see :ref:`installation`) or set `PYTHONPATH` to find the Lcapy source files.

- Then fire up your favourite python interpreter, for example, ipython:

  >>> ipython --pylab


Simple circuit elements
=======================

The basic circuit elements are two-terminal (one-port) devices:

- Vdc DC voltage source

- Idc DC current source

- Vac AC voltage source

- Iac AC current source

- R resistance

- G conductance

- C capacitance

- L inductance

These are augmented by generic s-domain components:

- Y admittance

- Z impedance

- V voltage source

- I current source


Here are some examples of their creation:

   >>> from lcapy import *
   >>> R1 = R(10)
   >>> C1 = C(10e-6)
   >>> L1 = L('L_1')

Here a symbolic inductance is created (using the Python string
notation).  In each case, a reference to the circuit element object is
stored in a Python variable.  These can be printed using `print` or
`pprint` (pretty print), for example,

   >>> print(R1)
   R(10)
   >>> print(C1)
   print(C)
   C(1.00000000000000e-6)
   >>> print(L1)
   L(L_1)
   >>> pprint(L1)
   L(L₁)



Simple circuit element combinations
-----------------------------------

Here's an example of resistors in series

   >>> from lcapy import *
   >>> R1 = R(10)
   >>> R2 = R(5)
   >>> Rtot = R1 + R2
   >>> print(Rtot)
   R(10) + R(5)
   >>> print(Rtot.simplify())
   R(15)

Here `R(10)` creates a 10 ohm resistor and this is assigned to the
variable `R1`.  Similarly, `R(5)` creates a 5 ohm resistor and this is
assigned to the variable `R2`.  `Rtot` is the name of the network
formed by connecting `R1` and `R2` in series.  Calling the `simplify`
method will simplify the network and combine the resistors into a
single resistor equivalent.


Here's an example of a parallel combination of resistors.  Note that
the parallel operator is `|` instead of the usual `||`.

   >>> from lcapy import *
   >>> Rtot = R(10) | R(5)
   >>> print(Rtot)
   R(10) | R(5)
   >>> print(Rtot.simplify())
   R(10/3)

The result can be performed symbolically, for example,

   >>> from lcapy import *
   >>> Rtot = R('R_1') | R('R_2')
   >>> print(Rtot)
   R(R_1) | R(R_2)
   >>> print(Rtot.simplify())
   R(R_1*R_2/(R_1 + R_2))
   >>> Rtot.simplify().pprint()
   R(R₁) | R(R₂)

Notice the difference between `print` and `pprint` (pretty print) methods.

Here's another example using inductors in series

   >>> from lcapy import *
   >>> L1 = L(10)
   >>> L2 = L(5)
   >>> Ltot = L1 + L2
   >>> Ltot.pprint()
   L(10) + L(5)
   >>> Ltot.simplify().pprint()
   L(15)


Finally, here's an example of a parallel combination of capacitors

   >>> from lcapy import *
   >>> Ctot = C(10) | C(5)
   >>> Ctot.pprint()
   C(10) | C(5)
   >>> Ctot.simplify().pprint()
   C(15)


Impedances
----------

Let's consider a series R-L-C network

   >>> from lcapy import *
   >>> N = R(5) + L(20) + C(10)
   >>> N.pprint()
   R(5) + L(20) + C(10)
   >>> N.Z.pprint()
        2           
   200⋅s  + 80⋅s + 1
   ─────────────────
          20⋅s    

Notice the result is a rational function of `s`.  Remember impedance
is a frequency domain concept.  A rational function can be formatted
in a number of different ways, for example,

   >>> N.Z.ZPK().pprint()
      ⎛      ____    ⎞ ⎛      ____    ⎞
      ⎜    ╲╱ 14    1⎟ ⎜    ╲╱ 14    1⎟
   10⋅⎜s - ────── + ─⎟⋅⎜s + ────── + ─⎟
      ⎝      20     5⎠ ⎝      20     5⎠
   ────────────────────────────────────
                    s                 
   >>> N.Z.canonical().pprint()
      ⎛ 2   2⋅s    1 ⎞
   10⋅⎜s  + ─── + ───⎟
      ⎝      5    200⎠
   ───────────────────
            s         

Here `ZPK()` prints the impedance in ZPK (zero-pole-gain) form while
`canonical()` prints the numerator and denominator of the rational
function in monic form (with unity leading coefficient).

The corresponding parallel R-L-C network yields

   >>> from lcapy import *
   >>> N = R(5) | L(20) | C(10)
   >>> N.pprint()
   R(5) | L(20) | C(10)
   >>> N.Z.pprint()
         20⋅s      
   ────────────────
        2          
   200⋅s  + 4⋅s + 1

   >>> N.Z.ZPK().pprint()
                   s                 
   ──────────────────────────────────
      ⎛     1    7⋅ⅈ⎞ ⎛     1    7⋅ⅈ⎞
   10⋅⎜s + ─── - ───⎟⋅⎜s + ─── + ───⎟
      ⎝    100   100⎠ ⎝    100   100⎠
   >>> N.Z.canonical().pprint()
           s         
   ──────────────────
      ⎛ 2   s     1 ⎞
   10⋅⎜s  + ── + ───⎟
      ⎝     50   200⎠
   >>> N.Y.pprint()
        2          
   200⋅s  + 4⋅s + 1
   ────────────────
         20⋅s      

Notice how `N.Y` returns the admittance of the network, the reciprocal
of the impedance `N.Z`.


The frequency response can be evaluated numerically by specifying a
vector of time values.

   >>> from lcapy import *
   >>> from numpy import logspace
   >>> N = Vdc(20) + R(5) + C(10)
   >>> f = linspace(0, 4, 400)
   >>> Isc = N.Isc.frequency_response(f)

Then the frequency response can be plotted.  For example,

   >>> from matplotlib.pyplot import figure, show
   >>> fig = figure()
   >>> ax = fig.add_subplot(111)
   >>> ax.loglog(f, abs(Isc), linewidth=2)
   >>> ax.set_xlabel('Frequency (Hz)')
   >>> ax.set_ylabel('Current (A/Hz)')
   >>> ax.grid(True)
   >>> show()

.. image:: examples/series-VRC1-Isc.png
   :width: 15cm


Here's a complete example Python script to plot the impedance of a
series R-L-C network:

.. literalinclude:: examples/series-RLC3-Z.py


.. image:: examples/series-RLC3-Z.png
   :width: 15cm


Simple transient analysis
=========================

Let's consider a series R-C network in series with a DC voltage source
(well, its not really a DC voltage source but a DC voltage multiplied
by a unit step)

   >>> from lcapy import *
   >>> N = Vdc(20) + R(5) + C(10)
   >>> N.pprint()
   Vdc(20) + R(5) + C(10)
   >>> Voc = N.Voc
   >>> Voc.pprint()
   20
   ──
   s 
   >>> N.Isc.pprint()
   200   
   ────────
   50⋅s + 1
   >>> isc = N.Isc.transient_response()
   >>> isc.pprint()
      -t              
      ───             
       50             
   4⋅ℯ   ⋅Heaviside(t)

Here `N` is network formed by the components in series, and `N.Voc` is
the open-circuit voltage across the network.  Note, this is the same
as the s-domain value of the voltage source.  `N.Isc` is the
short-circuit s-domain voltage through the network.  The method
`transient_response` converts this to the time-domain, a damped
exponential response multiplied by the ubiquitous unit step.

Of course, this could have been done symbolically,

   >>> from lcapy import *
   >>> N = Vdc('V_1') + R('R_1') + C('C_1')
   >>> N.pprint()
   Vdc(V₁) + R(R₁) + C(C₁)
   >>> Voc = N.Voc
   >>> Voc.pprint()
   V₁
   ──
   s 
   >>> N.Isc.pprint()
   C₁⋅V₁   
   ───────────
   C₁⋅R₁⋅s + 1
   >>> isc = N.Isc.transient_response()
   >>> isc.pprint()
         -t               
        ─────             
        C₁⋅R₁             
    V₁⋅ℯ     ⋅Heaviside(t)
    ──────────────────────
          R₁          

The transient response can be evaluated numerically by specifying a
vector of time values.

   >>> from lcapy import *
   >>> from numpy import linspace
   >>> N = Vdc(20) + R(5) + C(10)
   >>> t = linspace(0, 100, 400)
   >>> isc = N.Isc.transient_response(t)

Then the transient response can be plotted.  For example,

   >>> from matplotlib.pyplot import figure, show
   >>> fig = figure()
   >>> ax = fig.add_subplot(111)
   >>> ax.plot(t, isc, linewidth=2)
   >>> ax.set_xlabel('Time (s)')
   >>> ax.set_ylabel('Current (A)')
   >>> ax.grid(True)
   >>> show()

.. image:: examples/series-VRC1-isc.png
   :width: 15cm


Here's a complete example Python script of the short-circuit current
through an underdamped series RLC network:

.. literalinclude:: examples/series-VRLC1-isc.py

.. image:: examples/series-VRLC1-isc.png
   :width: 15cm



Transformations
===============

A one-port network can be represented as a Thevenin network (a series
combination of a voltage source and an impedance) or as a Norton
network (a parallel combination of a current source and an
admittance).

Here's an example of a Thevenin to Norton transformation:

   >>> from lcapy import *
   >>> T = Vdc(10) + R(5)
   >>> N = T.norton()
   >>> N.pprint()
   G(1/5) | Idc(2)

Similarly, here's an example of a Norton to Thevenin transformation:

   >>> from lcapy import *
   >>> N = Idc(10) | R(5)
   >>> T = N.thevenin()
   >>> T.pprint()
   R(5) + Vdc(50)



Two-port networks
=================

The basic circuit elements are one-port networks.  They can be
combined to create a two-port network.  The simplest two-port is a
shunt::

         -----+----
              |    
            +-+-+  
            |   |  
            |OP |  
            |   |  
            +-+-+  
              |    
         -----+----

A more interesting two-port network is an L section (voltage divider)::

           +---------+       
         --+   OP1   +---+----
           +---------+   |   
                       +-+-+ 
                       |   | 
                       |OP2| 
                       |   | 
                       +-+-+ 
                         |   
         ----------------+----

This is comprised from any two one-port networks.  For example,
   >>> from lcapy import *
   >>> R1 = R('R_1')
   >>> R2 = R('R_2')
   >>> N = LSection(R1, R2)
   >>> N.Vtransfer.pprint()
   R_2/(R_1 + R_2)

Here `N.Vtransfer` determines the forward voltage transfer function
`V_2(s) / V_1(s)`.

The open-circuit input impedance can be found using:
   >>> N.Z1oc.pprint()
   R₁ + R₂

The open-circuit output impedance can be found using:
   >>> N.Z2oc.pprint()
   R₂

The short-circuit input admittance can be found using:
   >>> N.Y1sc.pprint()
   1 
   ──
   R₁

The short-circuit output admittance can be found using:
   >>> N.Y2sc.pprint()
   R₁ + R₂
   ───────
    R₁⋅R₂ 


Two-port combinations
---------------------

Two-port networks can be combined in series, parallel, series at the
input with parallel at the output (hybrid), parallel at the input with
series at the output (inverse hybrid), but the most common is the
chain or cascade.  This connects the output of the first two-port to
the input of the second two-port.

For example, an L section can be created by chaining a shunt to a
series one-port.

   >>> from lcapy import *
   >>> N = Series(R('R_1')).chain(Shunt(R('R_2')))
   >>> N.Vtransfer.pprint()
   R_2/(R_1 + R_2)



Two-port matrices
-----------------

Two-port networks can be represented by six two by two matrices, A, B,
G, H, Y, Z.  Each has their own merits (see
http://en.wikipedia.org/wiki/Two-port_network).

Consider an L section comprised of two resistors:
   >>> from lcapy import *
   >>> N = LSection(R('R_1'), R('R_2')))

The different matrix representations can be shown using:
   >>> N.A.pprint()
   ⎡R₁ + R₂    ⎤
   ⎢───────  R₁⎥
   ⎢   R₂      ⎥
   ⎢           ⎥
   ⎢  1        ⎥
   ⎢  ──     1 ⎥
   ⎣  R₂       ⎦
   >>> N.B.pprint()
   ⎡ 1    -R₁  ⎤
   ⎢           ⎥
   ⎢-1   R₁    ⎥
   ⎢───  ── + 1⎥
   ⎣ R₂  R₂    ⎦
   >>> N.G.pprint()
   ⎡   1       -R₂  ⎤
   ⎢───────  ───────⎥
   ⎢R₁ + R₂  R₁ + R₂⎥
   ⎢                ⎥
   ⎢   R₂     R₁⋅R₂ ⎥
   ⎢───────  ───────⎥
   ⎣R₁ + R₂  R₁ + R₂⎦
   >>> N.H.pprint()
   ⎡R₁  1 ⎤
   ⎢      ⎥
   ⎢    1 ⎥
   ⎢-1  ──⎥
   ⎣    R₂⎦
   >>> N.Y.pprint()
   ⎡1      -1   ⎤
   ⎢──     ───  ⎥
   ⎢R₁      R₁  ⎥
   ⎢            ⎥
   ⎢-1   R₁ + R₂⎥
   ⎢───  ───────⎥
   ⎣ R₁   R₁⋅R₂ ⎦
   >>> N.Z.pprint()
   ⎡R₁ + R₂  R₂⎤
   ⎢           ⎥
   ⎣  R₂     R₂⎦


Note, some of the two-port matrices cannot represent a network.  For
example, a series impedance has a non specified Z matrix and a shunt
impedance has a non specified Y matrix.


Transfer Functions
==================

Transfer functions can be created in a similar manner to Matlab,
either using lists of numerator and denominator coefficients:

    >>> from lcapy import *
    >>>
    >>> H1 = tf(0.001, [1, 0.05, 0])
    >>> H1.pprint()
        0.001     
    ───────────────
         2         
    1.0⋅s  + 0.05⋅s

from lists of poles and zeros (and optional gain):

   >>> from lcapy import *
   >>>
   >>> H2 = zp2tf([], [0, -0.05])
   >>> H2.pprint()
        0.001     
   ───────────────
        2         
   1.0⋅s  + 0.05⋅s

or symbolically:

   >>> from lcapy import *
   >>>
   >>> H3 = 0.001 / (s**2 + 0.05 * s)
   >>> H3.pprint()
        0.001     
   ───────────────
        2         
   1.0⋅s  + 0.05⋅s


In each case, parameters can be expressed numerically or symbolically,
for example,

   >>> from lcapy import *
   >>>
   >>> H4 = zp2tf(['z_1'], ['p_1', 'p_2'])
   >>> H4.pprint()
          s - z₁      
   ───────────────────
   (-p₁ + s)⋅(-p₂ + s)



Partial Fraction Analysis
=========================

Lcapy can be used for converting rational functions into partial
fraction form.  Here's an example:

   >>> from lcapy import *
   >>>
   >>> G = 1 / (s**2 + 5 * s + 6)
   >>>
   >>> G.partfrac().pprint()
      1       1  
   - ───── + ─────
     s + 3   s + 2

Here's an example of a not strictly proper rational function,

   >>> from lcapy import *
   >>>
   >>> H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6)
   >>>
   >>> H.partfrac().pprint()
         70      90 
   5 + ───── - ─────
       s + 3   s + 2

The rational function can also be printed in ZPK form:

   >>> H.ZPK().pprint()
   5⋅(s - 4)⋅(s + 5)
   ─────────────────
    (s + 2)⋅(s + 3) 

Here it is obvious that the poles are -2 and -3.  These can also be found using the poles function:

   >>> H.poles()
   {-3: 1, -2: 1}

Here the number after the colon indicates how many times the pole is repeated.

Similarly, the zeros can be found using the zeros function:

   >>> H.zeros()
   {-5: 1, 4: 1}

Lcapy can also handle rational functions with a delay.



Inverse Laplace transforms
==========================

Lcapy can perform inverse Laplace transforms.   Here's an example for
a strictly proper rational function:

   >>> from lcapy import s
   >>> H = 5 * (s - 4) / (s**2 + 5 * s + 6)
   >>> H.partfrac().pprint()
     35      30 
   ───── - ─────
   s + 3   s + 2
   >>> H.inverse_laplace().pprint()
   ⎛      -2⋅t       -3⋅t⎞             
   ⎝- 30⋅ℯ     + 35⋅ℯ    ⎠⋅Heaviside(t)

The Heaviside function is the unit step.

When the rational function is not strictly proper, the inverse Laplace
transform has Dirac deltas (and derivatives of Dirac deltas):

   >>> from lcapy import s
   >>> H = 5 * (s - 4) / (s**2 + 5 * s + 6)
   >>> H.partfrac().pprint()
        70      90 
   5 + ───── - ─────
       s + 3   s + 2
   >>> H.inverse_laplace().pprint()
   ⎛      -2⋅t       -3⋅t⎞                               
   ⎝- 90⋅ℯ     + 70⋅ℯ    ⎠⋅Heaviside(t) + 5⋅DiracDelta(t)


Here's another example of a strictly proper rational function with a
repeated pole:

   >>> from lcapy import s
   >>> H = 5 * (s + 5) / ((s + 3) * (s + 3))
   >>> H.ZPK().pprint()
   5⋅(s + 5)
   ─────────
           2
    (s + 3) 
   >>> H.partfrac().pprint()
     5        10   
   ───── + ────────
   s + 3          2
           (s + 3) 
   >>> H.inverse_laplace().pprint()
   ⎛      -3⋅t      -3⋅t⎞             
   ⎝10⋅t⋅ℯ     + 5⋅ℯ    ⎠⋅Heaviside(t)


Rational functions with delays can also be handled:

   >>> from lcapy import s
   >>> import sympy as sym
   >>> T = sym.symbols('T')
   >>> H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6) * sym.exp(-s * T)
   >>> H.partfrac().pprint()
   ⎛      70      90 ⎞  -T⋅s
   ⎜5 + ───── - ─────⎟⋅ℯ    
   ⎝    s + 3   s + 2⎠      
   >>> H.inverse_laplace().pprint()
   ⎛      2⋅T - 2⋅t       3⋅T - 3⋅t⎞                                         
   ⎝- 90⋅ℯ          + 70⋅ℯ         ⎠⋅Heaviside(-T + t) + 5⋅DiracDelta(-T + t)



Laplace transforms
==================

Lcapy can also perform Laplace transforms.   Here's an example:

   >>> from lcapy import t
   >>> 
   >>> v = 10 * t ** 2 + 3 * t
   >>> v.laplace().pprint()
   3⋅s + 20
   ────────
       3   
      s   


Circuit analysis
================

The nodal voltages for a linear circuit can be found using Modified
Nodal Analysis (MNA).  This requires the circuit topology be entered as
a netlist.  This describes each component, its name, value, and the
nodes it is connected to.  This netlist can be read from a file or
created dynamically, for example,

   >>> from lcapy import Circuit
   >>> cct = Circuit()
   >>> cct.add('V1 1 0 dc 10') 
   >>> cct.add('Ra 1 2 3e3') 
   >>> cct.add('Rb 2 0 1e3') 

This creates a circuit comprised of a 10 V DC voltage source connected
to two resistors in series.  The node named 0 denotes the ground which
the other voltages are referenced to.

The circuit has a number of attributes that can be interrogated to
find circuit voltages and currents:

- `V` s-domain node voltage directory indexed by node name

- `I` s-domain branch current directory indexed by component name

- `Vd` s-domain branch voltage difference directory indexed by component name

- `v` t-domain node voltage directory indexed by node name

- `i` t-domain branch current directory indexed by component name

- `vd` t-domain branch voltage difference directory indexed by component name


For example,

   >>> cct.V[1]
   10.0
   ────
    s  
   >>> cct.V[2]
   2.5
   ───
    s 

Notice, how the displayed voltages are Laplace domain voltages.  The
transient voltages can be determined using an inverse Laplace transform:

   >>> cct.V[1].inverse_laplace().pprint()
   10.0⋅Heaviside(t)

Alternatively, 

   >>> cct.v[1]
   10.0⋅Heaviside(t)


For another example, the s-domain voltage difference across the
resistor Ra can be found using:

   >>> cct.Vd['Ra'].pprint()
   7.5
   ───
    s 

Since Lcapy uses SymPy, circuit analysis can be performed
symbolically.  This can be achieved by using symbolic arguments or by
not specifying a component value.  In the latter case, Lcapy will
use the component name for its value.  For example,

   >>> cct = Circuit()
   >>> cct.add('V1 1 0 dc Vs') 
   >>> cct.add('R1 1 2') 
   >>> cct.add('C1 2 0') 
   >>> cct.V[2].pprint()
        Vs     
   ────────────
          2    
   C₁⋅R₁⋅s  + s
   >>> : cct.V[2].inverse_laplace().pprint()
   ⎛          -t  ⎞             
   ⎜         ─────⎟             
   ⎜         C₁⋅R₁⎟             
   ⎝Vs - Vs⋅ℯ     ⎠⋅Heaviside(t)


Initial Conditions
------------------

The initial voltage difference across a capacitor or the initial
current through an inductor can be specified as the second argument.
For example,

   >>> cct = Circuit()
   >>> cct.add('V1 1 0 dc Vs') 
   >>> cct.add('C1 2 1 C1 v0') 
   >>> cct.add('L1 2 0 L1 i0') 
   >>> cct.V[2].pprint()
   C₁⋅L₁⋅Vs⋅s + C₁⋅L₁⋅s⋅v₀ - L₁⋅i₀
   ───────────────────────────────
                    2             
             C₁⋅L₁⋅s  + 1   


Transfer functions
------------------

Transfer functions can be found from the ratio of two s-domain
quantities such as voltage or current with zero initial conditions.
Here's an example using an arbitrary input voltage `V(s)`:

   >>> from lcapy import Circuit
   >>> cct = Circuit()
   >>> cct.add('V1 1 0 V(s)') 
   >>> cct.add('R1 1 2') 
   >>> cct.add('C1 2 0') 
   >>> cct.V[2].pprint()
       V(s)   
   ───────────
   C₁⋅R₁⋅s + 1

   >>> H = cct.V[2] / cct.V[1]
   >>> H.pprint()
        1     
   ───────────
   C₁⋅R₁⋅s + 1

The corresponding impulse response can found from an inverse Laplace transform:

   >>> H.inverse_laplace().pprint()
     -t               
    ─────             
    C₁⋅R₁             
   ℯ     ⋅Heaviside(t)
   ───────────────────
          C₁⋅R₁ 


Netlist specification
---------------------

The general form for a net is:

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

- Arbitrary time domain voltage source:

   vname Np Nm vexpr

- DC current source of current I:

   Iname Np Nm dc I

- AC current source of current I, frequency f, and phase p:

   Iname Np Nm ac I p

- Arbitrary s-domain current source:

   Iname Np Nm Iexpr

- Arbitrary time domain current source:

   iname Np Nm iexpr

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

Here's an example of a ramp voltage source (note, the expression cannot have
spaces in it):

v1 1 0 t*Heaviside(t)

Here's an example of a Dirac Delta voltage source:

v2 1 0 DiracDelta(t)

Here's an example of a cosine current current of amplitude 20 A and frequency f

i1 1 0 20*cos(2*pi*f*t)

Here's an example of a negative exponential current of amplitude 20 A
and time constant 2 s

i1 1 0 20*exp(-t/4)

Here's another way of specifying a 20 V DC voltage source, this time
using its Laplace transform:

V1 1 0 20/s



Other methods
-------------

   cct.Isc(n1, n2)      Short-circuit s-domain current between nodes n1 and n2.

   cct.Voc(n1, n2)      Open-circuit s-domain voltage between nodes n1 and n2.

   cct.isc(n1, n2)      Short-circuit t-domain current between nodes n1 and n2.

   cct.voc(n1, n2)      Open-circuit t-domain voltage between nodes n1 and n2.
   
   cct.Y(n1, n2)        Admittance between nodes n1 and n2.
  
   cct.Z(n1, n2)        Impedance between nodes n1 and n2.

   cct.kill()           Remove independent sources.

   cct.transfer(n1, n2, n3, n4) Transfer function for voltage difference between nodes n3 and n4 compared to the voltage difference between nodes n1 and n2.

   cct.thevenin(n1, n2) Thevenin model between nodes n1 and n2.

   cct.norton(n1, n2)    Norton model between nodes n1 and n2.

   cct.twoport(self, n1, n2, n3, n4) Create two-port component where
        I1 is the current flowing into n1 and out of n2, I2 is the
        current flowing into n3 and out of n4, V1 = V[n1] - V[n2], V2
        = V[n3] - V[n4].

   cct.remove(component) Remove component from net list.

   cct.netfile_add(filename) Add netlist from file.

   cct.s_model()         Convert circuit to s-domain model.

   cct.pre_initital_model()   Convert circuit to pre-initial model.

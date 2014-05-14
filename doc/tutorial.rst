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

- You will need to install SymPy, see http://docs.sympy.org/latest/install.html.

- For plotting you will need to install matplotlib.

- Lcapy can be downloaded from https://github.com/mph-/lcapy

  >>> git clone https://github.com/mph-/lcapy

- Lcapy can be installed using:

  >>> cd lcapy
  >>> sudo python setup.py install

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
===================================

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
   >>> pprint(Rtot.simplify())
   R(R₁) | R(R₂)

Notice the difference between `print` and `pprint` (pretty print).

Here's another example using inductors in series

   >>> from lcapy import *
   >>> L1 = L(10)
   >>> L2 = L(5)
   >>> Ltot = L1 + L2
   >>> pprint(Ltot)
   L(10) + L(5)
   >>> pprint(Ltot.simplify())
   L(15)


Finally, here's an example of a parallel combination of capacitors

   >>> from lcapy import *
   >>> Ctot = C(10) | C(5)
   >>> pprint(Ctot)
   C(10) | C(5)
   >>> pprint(Ctot.simplify())
   C(15)


Impedances
==========

Let's consider a series R-L-C network

   >>> from lcapy import *
   >>> N = R(5) + L(20) + C(10)
   >>> pprint(N)
   R(5) + L(20) + C(10)
   >>> pprint(N.Z)
        2           
   200⋅s  + 80⋅s + 1
   ─────────────────
          20⋅s    

Notice the result is a rational function of `s`.  Remember impedance
is a frequency domain concept.  A rational function can be formatted
in a number of different ways, for example,

   >>> pprint(N.Z.ZPK())
      ⎛      ____    ⎞ ⎛      ____    ⎞
      ⎜    ╲╱ 14    1⎟ ⎜    ╲╱ 14    1⎟
   10⋅⎜s - ────── + ─⎟⋅⎜s + ────── + ─⎟
      ⎝      20     5⎠ ⎝      20     5⎠
   ────────────────────────────────────
                    s                 
   >>> pprint(N.Z.canonical())
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
   >>> pprint(N)
   R(5) | L(20) | C(10)
   >>> pprint(N.Z)
         20⋅s      
   ────────────────
        2          
   200⋅s  + 4⋅s + 1

   >>> pprint(N.Z.ZPK())
                   s                 
   ──────────────────────────────────
      ⎛     1    7⋅ⅈ⎞ ⎛     1    7⋅ⅈ⎞
   10⋅⎜s + ─── - ───⎟⋅⎜s + ─── + ───⎟
      ⎝    100   100⎠ ⎝    100   100⎠
   >>> pprint(N.Z.canonical())
           s         
   ──────────────────
      ⎛ 2   s     1 ⎞
   10⋅⎜s  + ── + ───⎟
      ⎝     50   200⎠
   >>> pprint(N.Y)
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
   >>> pprint(N)
   Vdc(20) + R(5) + C(10)
   >>> Voc = N.Voc
   >>> pprint(Voc)
   20
   ──
   s 
   >>> pprint(N.Isc)
   200   
   ────────
   50⋅s + 1
   >>> isc = N.Isc.transient_response()
   >>> pprint(isc)
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
   >>> pprint(N)
   Vdc(V₁) + R(R₁) + C(C₁)
   >>> Voc = N.Voc
   >>> pprint(Voc)
   V₁
   ──
   s 
   >>> pprint(N.Isc)
   C₁⋅V₁   
   ───────────
   C₁⋅R₁⋅s + 1
   >>> isc = N.Isc.transient_response()
   >>> pprint(isc)
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
   >>> pprint(N)
   G(1/5) | Idc(2)

Similarly, here's an example of a Norton to Thevenin transformation:

   >>> from lcapy import *
   >>> N = Idc(10) | R(5)
   >>> T = N.thevenin()
   >>> pprint(T)
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
   >>> pprint(N.Vtransfer)
   R_2/(R_1 + R_2)

Here `N.Vtransfer` determines the forward voltage transfer function
`V_2(s) / V_1(s)`.

The open-circuit input impedance can be found using:
   >>> pprint(N.Z1oc)
   R₁ + R₂

The open-circuit output impedance can be found using:
   >>> pprint(N.Z2oc)
   R₂

The short-circuit input admittance can be found using:
   >>> pprint(N.Y1sc)
   1 
   ──
   R₁

The short-circuit output admittance can be found using:
   >>> pprint(N.Y2sc)
   R₁ + R₂
   ───────
    R₁⋅R₂ 


Two-port combinations
=====================

Two-port networks can be combined in series, parallel, series at the
input with parallel at the output (hybrid), parallel at the input with
series at the output (inverse hybrid), but the most common is the
chain or cascade.  This connects the output of the first two-port to
the input of the second two-port.

For example, an L section can be created by chaining a shunt to a
series one-port.

   >>> from lcapy import *
   >>> N = Series(R('R_1')).chain(Shunt(R('R_2')))
   >>> pprint(N.Vtransfer)
   R_2/(R_1 + R_2)



Two-port matrices
=================

Two-port networks can be represented by six two by two matrices, A, B,
G, H, Y, Z.  Each has their own merits (see
http://en.wikipedia.org/wiki/Two-port_network).

Consider an L section comprised of two resistors:
   >>> from lcapy import *
   >>> N = LSection(R('R_1'), R('R_2')))

The different matrix representations can be shown using:
   >>> pprint(N.A)
   ⎡R₁ + R₂    ⎤
   ⎢───────  R₁⎥
   ⎢   R₂      ⎥
   ⎢           ⎥
   ⎢  1        ⎥
   ⎢  ──     1 ⎥
   ⎣  R₂       ⎦
   >>> pprint(N.B)
   ⎡ 1    -R₁  ⎤
   ⎢           ⎥
   ⎢-1   R₁    ⎥
   ⎢───  ── + 1⎥
   ⎣ R₂  R₂    ⎦
   >>> pprint(N.G)
   ⎡   1       -R₂  ⎤
   ⎢───────  ───────⎥
   ⎢R₁ + R₂  R₁ + R₂⎥
   ⎢                ⎥
   ⎢   R₂     R₁⋅R₂ ⎥
   ⎢───────  ───────⎥
   ⎣R₁ + R₂  R₁ + R₂⎦
   >>> pprint(N.H)
   ⎡R₁  1 ⎤
   ⎢      ⎥
   ⎢    1 ⎥
   ⎢-1  ──⎥
   ⎣    R₂⎦
   >>> pprint(N.Y)
   ⎡1      -1   ⎤
   ⎢──     ───  ⎥
   ⎢R₁      R₁  ⎥
   ⎢            ⎥
   ⎢-1   R₁ + R₂⎥
   ⎢───  ───────⎥
   ⎣ R₁   R₁⋅R₂ ⎦
   >>>pprint(N.Z)
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
    >>> pprint(H1)
        0.001     
    ───────────────
         2         
    1.0⋅s  + 0.05⋅s

from lists of poles and zeros (and optional gain):

   >>> from lcapy import *
   >>>
   >>> H2 = zp2tf([], [0, -0.05])
   >>> pprint(H2)
        0.001     
   ───────────────
        2         
   1.0⋅s  + 0.05⋅s

or symbolically:

   >>> from lcapy import *
   >>>
   >>> s = tf('s')
   >>> H3 = 0.001 / (s**2 + 0.05 * s)
   >>> pprint(H3)
        0.001     
   ───────────────
        2         
   1.0⋅s  + 0.05⋅s


In each case, parameters can be expressed numerically or symbolically, for example,

   >>> from lcapy import *
   >>>
   >>> H4 = zp2tf(['z_1'], ['p_1', 'p_2'])
   >>> pprint(H4)
          s - z₁      
   ───────────────────
   (-p₁ + s)⋅(-p₂ + s)



Partial Fraction Analysis
=========================

Lcapy can be used for converting rational functions into partial
fraction form.  Here's an example:

   >>> from lcapy import *
   >>>
   >>> s = sExpr.s
   >>>
   >>> G = 1 / (s**2 + 5 * s + 6)
   >>>
   >>> pprint(partfrac(G))
      1       1  
   - ───── + ─────
     s + 3   s + 2

Here's an example of a not strictly proper rational function,

   >>> from lcapy import *
   >>>
   >>> s = sExpr.s
   >>>
   >>> H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6)
   >>>
   >>> pprint(partfrac(H))
         70      90 
   5 + ───── - ─────
       s + 3   s + 2

The rational function can also be printed in ZPK form:

   >>> pprint(ZPK(H))
   5⋅(s - 4)⋅(s + 5)
   ─────────────────
    (s + 2)⋅(s + 3) 

Here it is obvious that the poles are -2 and -3.  These can also be found using the poles function:

   >>> poles(H)
   {-3: 1, -2: 1}

Here the number after the colon indicates how many times the pole is repeated.

Similarly, the zeros can be found using the zeros function:

   >>> zeros(H)
   {-5: 1, 4: 1}

Lcapy can also handle rational functions with a delay.



Inverse Laplace transforms
==========================

Lcapy can perform inverse Laplace transforms.   Here's an example for
a strictly proper rational function:

   >>> from lcapy import *
   >>> s = sExpr.s
   >>> H = 5 * (s - 4) / (s**2 + 5 * s + 6)
   >>> pprint(partfrac(H))
     35      30 
   ───── - ─────
   s + 3   s + 2
   >>> pprint(inverse_laplace(H))
   ⎛      -2⋅t       -3⋅t⎞             
   ⎝- 30⋅ℯ     + 35⋅ℯ    ⎠⋅Heaviside(t)

The Heaviside function is the unit step.

When the rational function is not strictly proper, the inverse Laplace
transform has Dirac deltas (and derivatives of Dirac deltas):

   >>> from lcapy import *
   >>> s = sExpr.s
   >>> H = 5 * (s - 4) / (s**2 + 5 * s + 6)
   >>> pprint(partfrac(H)) 
        70      90 
   5 + ───── - ─────
       s + 3   s + 2
   >>> pprint(inverse_laplace(H))
   ⎛      -2⋅t       -3⋅t⎞                               
   ⎝- 90⋅ℯ     + 70⋅ℯ    ⎠⋅Heaviside(t) + 5⋅DiracDelta(t)


Here's another example of a strictly proper rational function with a
repeated pole:

   >>> from lcapy import *
   >>> s = sExpr.s
   >>> H = 5 * (s + 5) / ((s + 3) * (s + 3))
   >>> pprint(ZPK(H))
   5⋅(s + 5)
   ─────────
           2
    (s + 3) 
   >>> pprint(partfrac(H))
     5        10   
   ───── + ────────
   s + 3          2
           (s + 3) 
   >>> pprint(inverse_laplace(H))
   ⎛      -3⋅t      -3⋅t⎞             
   ⎝10⋅t⋅ℯ     + 5⋅ℯ    ⎠⋅Heaviside(t)


Rational functions with delays can also be handled:

   >>> from lcapy import *
   >>> from sympy import symbols
   >>> s, T = sym.symbols('s T')
   >>> H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6) * sym.exp(-s * T)
   >>> pprint(partfrac(H))
   ⎛      70      90 ⎞  -T⋅s
   ⎜5 + ───── - ─────⎟⋅ℯ    
   ⎝    s + 3   s + 2⎠      
   >>> pprint(inverse_laplace(H))
   ⎛      2⋅T - 2⋅t       3⋅T - 3⋅t⎞                                         
   ⎝- 90⋅ℯ          + 70⋅ℯ         ⎠⋅Heaviside(-T + t) + 5⋅DiracDelta(-T + t)



Circuit analysis
================

The nodal voltages for a linear circuit can be found using Modified
Nodal Analysis (MNA).  This requires the circuit topology be entered as
a netlist.  This describes each component, its name, value, and the
nodes it is connected to.  This netlist can be read from a file or
created dynamically, for example

   >>> from lcapy import pprint, Circuit
   >>> cct = Circuit('Voltage divider')
   >>> cct.net_add('Vs 1 0 10') 
   >>> cct.net_add('Ra 1 2 3e3') 
   >>> cct.net_add('Rb 2 0 1e3') 

This creates a circuit comprised of a 10\,V DC voltage source connected
to two resistors in series.  The node named 0 denotes the ground which
the other voltages are referenced to.

The node voltages are stored in a directory indexed by the node name.  For example,
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

   >>> pprint(cct.V[1].inverse_laplace())
   10.0⋅Heaviside(t)


The branch voltages are stored in a directory indexed by component
name.  For example, the current through the voltage source Vs is:

   >>> pprint(cct.I['Vs'])
   0.1
   ───
    s  

Since Lcapy uses SymPy, circuit analysis can be performed
symbolically.  This can be achieved by not specifying a component
value.  Lcapy, will then create a symbol using the component name.

   >>> cct = Circuit('Series V R C')
   >>> cct.net_add('Vs 1 0') 
   >>> cct.net_add('R1 1 2') 
   >>> cct.net_add('C1 2 0') 
   >>> pprint(cct.V[2])
        Vs     
   ────────────
          2    
   C₁⋅R₁⋅s  + s
   >>> : pprint(cct.V[2].inverse_laplace())
   ⎛          -t  ⎞             
   ⎜         ─────⎟             
   ⎜         C₁⋅R₁⎟             
   ⎝Vs - Vs⋅ℯ     ⎠⋅Heaviside(t)

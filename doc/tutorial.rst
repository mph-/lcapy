.. TODO, split into manual and tutorial

========
Tutorial
========


Introduction
============

Lcapy is a Python package for linear circuit analysis.  It will only
solve linear, time invariant networks.  In other words, networks
comprised of basic circuit components (R, L, C, etc.) that do not vary
with time.

Lcapy can semi-automate the drawing of schematics from a netlist.

Lcapy cannot directly analyse non-linear devices such as diodes or
transistors although it does support simple opamps without saturation.
Nevertheless, it can draw them!

Lcapy uses Sympy (symbolic Python) for its values and expressions
and thus the circuit analysis can be performed symbolically.  See http://docs.sympy.org/latest/tutorial/index.html for the Sympy tutorial.

Internally, the circuit components are stored using their s-domain
equivalents, such as impedances and admittances.  This is convenient
for frequency response analysis but requires an inverse Laplace
transform for transient response analysis.


Preliminaries
=============

- Before you can use Lcapy you need to install the Lcapy package (see :ref:`installation`) or set `PYTHONPATH` to find the Lcapy source files.

- Then fire up your favourite python interpreter, for example, ipython:

  >>> ipython --pylab


Expressions
===========

Lcapy defines a number of symbols corresponding to different domains:

- t -- time (real)

- f -- frequency (real)

- s -- complex (s-domain) frequency

- omega -- angular frequency (real)

Expressions can be formed using these symbols, for example, a
time-domain expression can be created using:

   >>> from lcapy import t, delta, u
   >>> v = 2 * t * u(t) + 3 + delta(t)
   >>> i = 0 * t + 3

and a s-domain expression can be created using:

   >>> from lcapy import s, j, omega
   >>> H = (s + 3) / (s - 4)
   >>> H
   s + 3
   ─────
   s - 4

For steady-state signals, the s-domain can be converted to the angular
frequency domain by substituting :math:`\mathrm{j} \omega` for :math:`s`

   >>> from lcapy import s, j, omega
   >>> H = (s + 3) / (s - 4)
   >>> A = H(j * omega)
   >>> A
   ⅈ⋅ω + 3
   ───────
   ⅈ⋅ω - 4

Also note, real numbers are approximated by rationals.

Lcapy expressions have a number of attributes, including:

- numerator, N --  numerator of rational function

- denominator, D --  denominator of rational function

- magnitude -- magnitude

- angle -- angle

- real -- real part

- imag -- imaginary part

- conjugate -- complex conjugate

- expr -- the underlying Sympy expression

- val -- the expression as evaluated as a floating point value (if possible)

and a number of generic methods including:

- symplify -- attempt simple simplification of the expression

- rationalize_denominator -- multiply numerator and denominator by complex conjugate of denominator

- evaluate -- evaluate at specified vector and return floating point vector

Here's an example of using these attributes and methods:

   >>> from lcapy import s, j, omega
   >>> H = (s + 3) / (s - 4)
   >>> A = H(j * omega)
   >>> A
   ⅈ⋅ω + 3
   ───────
   ⅈ⋅ω - 4
   >>> A.rationalize_denominator()
    2             
   ω  - 7⋅ⅈ⋅ω - 12
   ───────────────
        2         
       ω  + 16  
   >>> A.real
    2     
   ω  - 12
   ───────
    2     
   ω  + 16
   >>> A.imag
    -7⋅ω  
   ───────
    2     
   ω  + 16
   >>> A.N
   ⅈ⋅ω + 3
   >>> A.D
   ⅈ⋅ω - 4
   >>> A.phase
        ⎛       2     ⎞
   atan2⎝-7⋅ω, ω  - 12⎠
   >>> A.magnitude
      __________________
     ╱  4       2       
   ╲╱  ω  + 25⋅ω  + 144 
   ─────────────────────
           2            
          ω  + 16       



Each domain has specific methods, including:

- fourier  

- laplace

- inverse_fourier

- inverse_laplace


Lcapy defines a number of functions that can be used in expressions, including:

- u --  Heaviside's unit step

- H -- Heaviside's unit step

- delta -- Dirac delta

- cos -- cosine

- sin -- sine

- sqrt -- square root

- exp -- exponential

- log10 -- logarithm base 10

- log -- natural logarithm




Simple circuit components
=========================

The basic circuit components are two-terminal (one-port) devices:

- I current source

- V voltage source

- R resistance

- C capacitance

- L inductance

These are augmented by generic s-domain components:

- Y admittance

- Z impedance


Here are some examples of their creation:

   >>> from lcapy import *
   >>> R1 = R(10)
   >>> C1 = C(10e-6)
   >>> L1 = L('L_1')



Simple circuit element combinations
-----------------------------------

Here's an example of resistors in series

   >>> from lcapy import *
   >>> R1 = R(10)
   >>> R2 = R(5)
   >>> Rtot = R1 + R2
   >>> Rtot
   R(10) + R(5)
   >>> Rtot.simplify()
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
   >>> Rtot
   R(10) | R(5)
   >>> Rtot.simplify()
   R(10/3)

The result can be performed symbolically, for example,

   >>> from lcapy import *
   >>> Rtot = R('R_1') | R('R_2')
   >>> Rtot
   R(R_1) | R(R_2)
   >>> Rtot.simplify()
   R(R_1*R_2/(R_1 + R_2))
   >>> Rtot.simplify()
   R(R₁) | R(R₂)

Here's another example using inductors in series

   >>> from lcapy import *
   >>> L1 = L(10)
   >>> L2 = L(5)
   >>> Ltot = L1 + L2
   >>> Ltot
   L(10) + L(5)
   >>> Ltot.simplify()
   L(15)


Finally, here's an example of a parallel combination of capacitors

   >>> from lcapy import *
   >>> Ctot = C(10) | C(5)
   >>> Ctot
   C(10) | C(5)
   >>> Ctot.simplify()
   C(15)


Impedances
----------

Let's consider a series R-L-C network

   >>> from lcapy import *
   >>> N = R(5) + L(20) + C(10)
   >>> N
   R(5) + L(20) + C(10)
   >>> N.Z
        2           
   200⋅s  + 80⋅s + 1
   ─────────────────
          20⋅s    

Notice the result is a rational function of `s`.  Remember impedance
is a frequency domain concept.  A rational function can be formatted
in a number of different ways, for example,

   >>> N.Z.ZPK()
      ⎛      ____    ⎞ ⎛      ____    ⎞
      ⎜    ╲╱ 14    1⎟ ⎜    ╲╱ 14    1⎟
   10⋅⎜s - ────── + ─⎟⋅⎜s + ────── + ─⎟
      ⎝      20     5⎠ ⎝      20     5⎠
   ────────────────────────────────────
                    s                 
   >>> N.Z.canonical()
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
   >>> N
   R(5) | L(20) | C(10)
   >>> N.Z
         20⋅s      
   ────────────────
        2          
   200⋅s  + 4⋅s + 1

   >>> N.Z.ZPK()
                   s                 
   ──────────────────────────────────
      ⎛     1    7⋅ⅈ⎞ ⎛     1    7⋅ⅈ⎞
   10⋅⎜s + ─── - ───⎟⋅⎜s + ─── + ───⎟
      ⎝    100   100⎠ ⎝    100   100⎠
   >>> N.Z.canonical()
           s         
   ──────────────────
      ⎛ 2   s     1 ⎞
   10⋅⎜s  + ── + ───⎟
      ⎝     50   200⎠
   >>> N.Y
        2          
   200⋅s  + 4⋅s + 1
   ────────────────
         20⋅s      

Notice how `N.Y` returns the admittance of the network, the reciprocal
of the impedance `N.Z`.


The frequency response can be evaluated numerically by specifying a
vector of frequency values.

   >>> from lcapy import *
   >>> from numpy import linspace
   >>> N = Vdc(20) + R(5) + C(10)
   >>> vf = linspace(0, 4, 400)
   >>> Isc = N.Isc.frequency_response().evaluate(vf)

Then the frequency response can be plotted.  For example,

   >>> from matplotlib.pyplot import figure, show
   >>> fig = figure()
   >>> ax = fig.add_subplot(111)
   >>> ax.loglog(f, abs(Isc), linewidth=2)
   >>> ax.set_xlabel('Frequency (Hz)')
   >>> ax.set_ylabel('Current (A/Hz)')
   >>> ax.grid(True)
   >>> show()

A simpler approach is to use the plot method:

   >>> from lcapy import *
   >>> from numpy import linspace
   >>> N = Vdc(20) + R(5) + C(10)
   >>> vf = linspace(0, 4, 400)
   >>> Isc = N.Isc.frequency_response().plot(vf, log_scale=True)

.. image:: examples/netlists/series-VRC1-Isc.png
   :width: 15cm

Here's a complete example Python script to plot the impedance of a
series R-L-C network:

.. literalinclude:: examples/netlists/series-RLC3-Z.py


.. image:: examples/netlists/series-RLC3-Z.png
   :width: 15cm


Simple transient analysis
=========================

Let's consider a series R-C network in series with a DC voltage source
(well, its not really a DC voltage source but a DC voltage multiplied
by a unit step)

   >>> from lcapy import *
   >>> N = Vdc(20) + R(5) + C(10)
   >>> N
   Vdc(20) + R(5) + C(10)
   >>> Voc = N.Voc
   >>> Voc
   20
   ──
   s 
   >>> N.Isc
   200   
   ────────
   50⋅s + 1
   >>> isc = N.Isc.transient_response()
   >>> isc
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
   >>> N
   Vdc(V₁) + R(R₁) + C(C₁)
   >>> Voc = N.Voc
   >>> Voc
   V₁
   ──
   s 
   >>> N.Isc
   C₁⋅V₁   
   ───────────
   C₁⋅R₁⋅s + 1
   >>> isc = N.Isc.transient_response()
   >>> isc
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

Then the transient response can be plotted.  Alternatively, the plot
method can be used.

.. literalinclude:: examples/netlists/series-VRC1-isc.py

This produces:

.. image:: examples/netlists/series-VRC1-isc.png
   :width: 15cm


Here's a complete example Python script of the short-circuit current
through an underdamped series RLC network:

.. literalinclude:: examples/netlists/series-VRLC1-isc.py

.. image:: examples/netlists/series-VRLC1-isc.png
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
   >>> N
   G(1/5) | Idc(2)

Similarly, here's an example of a Norton to Thevenin transformation:

   >>> from lcapy import *
   >>> N = Idc(10) | R(5)
   >>> T = N.thevenin()
   >>> T
   R(5) + Vdc(50)



Two-port networks
=================

The basic circuit components are one-port networks.  They can be
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
   >>> N.Vtransfer
   R_2/(R_1 + R_2)

Here `N.Vtransfer` determines the forward voltage transfer function
`V_2(s) / V_1(s)`.

The open-circuit input impedance can be found using:
   >>> N.Z1oc
   R₁ + R₂

The open-circuit output impedance can be found using:
   >>> N.Z2oc
   R₂

The short-circuit input admittance can be found using:
   >>> N.Y1sc
   1 
   ──
   R₁

The short-circuit output admittance can be found using:
   >>> N.Y2sc
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
   >>> N.Vtransfer
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
   >>> N.A
   ⎡R₁ + R₂    ⎤
   ⎢───────  R₁⎥
   ⎢   R₂      ⎥
   ⎢           ⎥
   ⎢  1        ⎥
   ⎢  ──     1 ⎥
   ⎣  R₂       ⎦
   >>> N.B
   ⎡ 1    -R₁  ⎤
   ⎢           ⎥
   ⎢-1   R₁    ⎥
   ⎢───  ── + 1⎥
   ⎣ R₂  R₂    ⎦
   >>> N.G
   ⎡   1       -R₂  ⎤
   ⎢───────  ───────⎥
   ⎢R₁ + R₂  R₁ + R₂⎥
   ⎢                ⎥
   ⎢   R₂     R₁⋅R₂ ⎥
   ⎢───────  ───────⎥
   ⎣R₁ + R₂  R₁ + R₂⎦
   >>> N.H
   ⎡R₁  1 ⎤
   ⎢      ⎥
   ⎢    1 ⎥
   ⎢-1  ──⎥
   ⎣    R₂⎦
   >>> N.Y
   ⎡1      -1   ⎤
   ⎢──     ───  ⎥
   ⎢R₁      R₁  ⎥
   ⎢            ⎥
   ⎢-1   R₁ + R₂⎥
   ⎢───  ───────⎥
   ⎣ R₁   R₁⋅R₂ ⎦
   >>> N.Z
   ⎡R₁ + R₂  R₂⎤
   ⎢           ⎥
   ⎣  R₂     R₂⎦


Note, some of the two-port matrices cannot represent a network.  For
example, a series impedance has a non specified Z matrix and a shunt
impedance has a non specified Y matrix.


Transfer functions
==================

Transfer functions can be created in a similar manner to Matlab,
either using lists of numerator and denominator coefficients:

    >>> from lcapy import *
    >>> H1 = tf(0.001, [1, 0.05, 0])
    >>> H1
        0.001     
    ───────────────
         2         
    1.0⋅s  + 0.05⋅s

from lists of poles and zeros (and optional gain):

   >>> from lcapy import *
   >>> H2 = zp2tf([], [0, -0.05])
   >>> H2
        0.001     
   ───────────────
        2         
   1.0⋅s  + 0.05⋅s

or symbolically:

   >>> from lcapy import *
   >>> H3 = 0.001 / (s**2 + 0.05 * s)
   >>> H3
        0.001     
   ───────────────
        2         
   1.0⋅s  + 0.05⋅s


In each case, parameters can be expressed numerically or symbolically,
for example,

   >>> from lcapy import *
   >>> H4 = zp2tf(['z_1'], ['p_1', 'p_2'])
   >>> H4
          s - z₁      
   ───────────────────
   (-p₁ + s)⋅(-p₂ + s)



Partial fraction analysis
=========================

Lcapy can be used for converting rational functions into partial
fraction form.  Here's an example:

   >>> from lcapy import *
   >>> G = 1 / (s**2 + 5 * s + 6)
   >>> G.partfrac()
      1       1  
   - ───── + ─────
     s + 3   s + 2

Here's an example of a not strictly proper rational function,

   >>> from lcapy import *
   >>> H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6)
   >>> H.partfrac()
         70      90 
   5 + ───── - ─────
       s + 3   s + 2

The rational function can also be printed in ZPK form:

   >>> H.ZPK()
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
   >>> H.partfrac()
     35      30 
   ───── - ─────
   s + 3   s + 2

   >>> H.inverse_laplace()
   ⎧      -2⋅t       -3⋅t           
   ⎨- 30⋅ℯ     + 35⋅ℯ      for t ≥ 0
   ⎩                                

Note that the unilateral inverse Laplace transform can only determine
the result for t >= 0.  If you know that the system is causal, then use


   >>> H.inverse_laplace(causal=True)

This gives
   ⎛      -2⋅t       -3⋅t⎞             
   ⎝- 30⋅ℯ     + 35⋅ℯ    ⎠⋅Heaviside(t)

The Heaviside function is the unit step.

When the rational function is not strictly proper, the inverse Laplace
transform has Dirac deltas (and derivatives of Dirac deltas):

   >>> from lcapy import s
   >>> H = 5 * (s - 4) / (s**2 + 5 * s + 6)
   >>> H.partfrac()
        70      90 
   5 + ───── - ─────
       s + 3   s + 2
   >>> H.inverse_laplace(causal=True)
   ⎛      -2⋅t       -3⋅t⎞                               
   ⎝- 90⋅ℯ     + 70⋅ℯ    ⎠⋅Heaviside(t) + 5⋅DiracDelta(t)


Here's another example of a strictly proper rational function with a
repeated pole:

   >>> from lcapy import s
   >>> H = 5 * (s + 5) / ((s + 3) * (s + 3))
   >>> H.ZPK()
   5⋅(s + 5)
   ─────────
           2
    (s + 3) 
   >>> H.partfrac()
     5        10   
   ───── + ────────
   s + 3          2
           (s + 3) 
   >>> H.inverse_laplace(causal=True)
   ⎛      -3⋅t      -3⋅t⎞             
   ⎝10⋅t⋅ℯ     + 5⋅ℯ    ⎠⋅Heaviside(t)


Rational functions with delays can also be handled:

   >>> from lcapy import s
   >>> import sympy as sym
   >>> T = sym.symbols('T')
   >>> H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6) * sym.exp(-s * T)
   >>> H.partfrac()
   ⎛      70      90 ⎞  -T⋅s
   ⎜5 + ───── - ─────⎟⋅ℯ    
   ⎝    s + 3   s + 2⎠      
   >>> H.inverse_laplace(causal=True)
   ⎛      2⋅T - 2⋅t       3⋅T - 3⋅t⎞                                         
   ⎝- 90⋅ℯ          + 70⋅ℯ         ⎠⋅Heaviside(-T + t) + 5⋅DiracDelta(-T + t)



Laplace transforms
==================

Lcapy can also perform Laplace transforms.   Here's an example:

   >>> from lcapy import t
   >>> v = 10 * t ** 2 + 3 * t
   >>> v.laplace()
   3⋅s + 20
   ────────
       3   
      s   


Circuit analysis
================

The nodal voltages for a linear circuit can be found using Modified
Nodal Analysis (MNA).  This requires the circuit topology be entered
as a netlist (see :ref:`netlists`).  This describes each component, its
name, value, and the nodes it is connected to.  This netlist can be
read from a file or created dynamically, for example,

   >>> from lcapy import Circuit
   >>> cct = Circuit()
   >>> cct.add('V1 1 0 dc 10') 
   >>> cct.add('Ra 1 2 3e3') 
   >>> cct.add('Rb 2 0 1e3') 

This creates a circuit comprised of a 10 V DC voltage source connected
to two resistors in series.  The node named 0 denotes the ground which
the other voltages are referenced to.

The circuit has an attribute for each circuit element (and for each
node starting with an alphabetical character).  These can be
interrogated to find the voltage drop across an element or the current
through an element, for example,

   >>> cct.V1.V
   10.0
   ────
    s 
   >>> cct.Rb.V
   2.5
   ───
    s  

Notice, how the displayed voltages are Laplace domain voltages.  The
transient voltages can be determined using an inverse Laplace transform:

   >>> cct.V1.V.inverse_laplace(causal=True)
   10.0⋅Heaviside(t)

Alternatively, using the lowercase `v` attribute:

   >>> cct.V1.v
   10.0⋅Heaviside(t)


The voltage between a node and ground can be determined with the node
name as an index, for example,

   >>> cct[1].V
   10.0
   ────
    s  
   >>> cct[2].V
   2.5
   ───
    s 

The circuit has a number of attributes that can be interrogated to
find circuit voltages and currents, including:

- `V` s-domain voltage directory indexed by node name or branch name

- `I` s-domain branch current directory indexed by component name

- `v` t-domain voltage directory indexed by node name or branch name

- `i` t-domain branch current directory indexed by component name


The circuit also has a number of methods, including:

- `Y` admittance between pair of nodes

- `Z` impedance between pair of nodes

- `Isc` short-circuit current between pair of nodes

- `Voc` open-circuit current between pair of nodes


Since Lcapy uses Sympy, circuit analysis can be performed
symbolically.  This can be achieved by using symbolic arguments or by
not specifying a component value.  In the latter case, Lcapy will
use the component name for its value.  For example,

   >>> cct = Circuit()
   >>> cct.add('V1 1 0 dc Vs') 
   >>> cct.add('R1 1 2') 
   >>> cct.add('C1 2 0') 
   >>> cct[2].V
        Vs     
   ────────────
          2    
   C₁⋅R₁⋅s  + s
   >>> : cct[2].V.inverse_laplace()
   ⎧            -t             
   ⎪           ─────           
   ⎨           C₁⋅R₁           
   ⎪V_s - V_s⋅ℯ       for t ≥ 0
   ⎩                           


Initial Conditions
------------------

The initial voltage difference across a capacitor or the initial
current through an inductor can be specified as the second argument.
For example,

   >>> cct = Circuit()
   >>> cct.add('V1 1 0 dc Vs') 
   >>> cct.add('C1 2 1 C1 v0') 
   >>> cct.add('L1 2 0 L1 i0') 
   >>> cct[2].V
   C₁⋅L₁⋅Vs⋅s + C₁⋅L₁⋅s⋅v₀ - L₁⋅i₀
   ───────────────────────────────
                    2             
             C₁⋅L₁⋅s  + 1   


Transfer functions
------------------

Transfer functions can be found from the ratio of two s-domain
quantities such as voltage or current with zero initial conditions.
Here's an example using an arbitrary input voltage `V(s)`

   >>> from lcapy import Circuit
   >>> cct = Circuit()
   >>> cct.add('V1 1 0 {V(s)}') 
   >>> cct.add('R1 1 2') 
   >>> cct.add('C1 2 0') 
   >>> cct[2].V
       V(s)   
   ───────────
   C₁⋅R₁⋅s + 1

   >>> H = cct[2].V / cct[1].V
   >>> H
        1     
   ───────────
   C₁⋅R₁⋅s + 1

The corresponding impulse response can found from an inverse Laplace transform:

   >>> H.inverse_laplace(causal=True)
     -t               
    ─────             
    C₁⋅R₁             
   ℯ     ⋅Heaviside(t)
   ───────────────────
          C₁⋅R₁ 

Transfer functions can also be created using the transfer method of a
circuit.  For example,

   >>> from lcapy import Circuit
   >>> cct = Circuit()
   >>> cct.add('R1 1 2') 
   >>> cct.add('C1 2 0') 
   >>> H = cct.transfer(1, 0, 2, 0)
   >>> H
        1     
   ───────────
   C₁⋅R₁⋅s + 1

In this example, the transfer method computes (V[1] - V[0]) / (V[2] -
V[0]).  In general, all independent sources are killed and so the
response is causal.

   >>> H.inverse_laplace()
     -t               
    ─────             
    C₁⋅R₁             
   ℯ     ⋅Heaviside(t)
   ───────────────────
          C₁⋅R₁       


Other circuit methods
---------------------

   cct.Isc(Np, Nm)      Short-circuit s-domain current between nodes Np and Nm.

   cct.Voc(Np, Nm)      Open-circuit s-domain voltage between nodes Np and Nm.

   cct.isc(Np, Nm)      Short-circuit t-domain current between nodes Np and Nm.

   cct.voc(Np, Nm)      Open-circuit t-domain voltage between nodes Np and Nm.
   
   cct.admittance(Np, Nm)        Admittance between nodes Np and Nm.
  
   cct.impedance(Np, Nm)         Impedance between nodes Np and Nm.

   cct.kill()           Remove independent sources.

   cct.kill_except(sources)      Remove independent sources except ones specified.

   cct.transfer(N1p, N1m, N2p, N2m) Voltage transfer function V2/V1, where V1 = V[N1p] - V[N1m], V2 = V[N2p] - V[N2m].

   cct.thevenin(Np, Nm) Thevenin model between nodes Np and Nm.

   cct.norton(Np, Nm)    Norton model between nodes Np and Nm.

   cct.twoport(self, N1p, N1m, N2p, N2m) Create two-port component where
        I1 is the current flowing into N1p and out of N1m, I2 is the
        current flowing into N2p and out of N2m, V1 = V[N1p] - V[N1m], V2
        = V[N2p] - V[N2m].

   cct.remove(component) Remove component from net list.

   cct.netfile_add(filename) Add netlist from file.

   cct.s_model()         Convert circuit to s-domain model.

   cct.pre_initial_model()   Convert circuit to pre-initial model.


Plotting
========

Lcapy expressions have a plot method; this differs depending on the
domain.  For example, the plot method for s-domain expressions
produces a pole-zero plot.  Here's an example:

.. literalinclude:: examples/netlists/tf1-pole-zero-plot.py

.. image:: examples/netlists/tf1-pole-zero-plot.png
   :width: 15cm


The plot method for f-domain and :math:`\omega` -domain expressions
produce spectral plots, for example, 


.. literalinclude:: examples/netlists/tf1-bode-plot.py

.. image:: examples/netlists/tf1-bode-plot.png
   :width: 15cm


Schematics
==========

Schematics can be generated from a netlist and from one port networks.
In both cases the drawing is performed using the LaTeX Circuitikz
package.  The schematic can be displayed interactively or saved to a
pdf or png file.


Netlist schematics
------------------

Hints are required to designate component orientation and explicit
wires are required to link nodes of the same potential but with
different coordinates.  For more details see :ref:`schematics`.

Here's an example:
   >>> from lcapy import Circuit
   >>> cct = Circuit()
   >>> cct.add('V1 1 0 {V(s)}; down') 
   >>> cct.add('R1 1 2; right') 
   >>> cct.add('C1 2 0_2; down') 
   >>> cct.add('W1 0 0_2; right') 
   >>> cct.draw('schematic.pdf')

Note, the orientation hints are appended to the netlist strings with a
semicolon delimiter.  The drawing direction is with respect to the
first node.  The component W1 is a wire.  Nodes with an underscore in
their name are not drawn with a closed blob.

Here's another example, this time loading the netlist from a file:
   >>> from lcapy import Circuit
   >>> cct = Circuit('voltage-divider.sch')
   >>> cct.draw('voltage-divider.pdf')

Here are the contents of the file 'voltage-divider.sch'::

   V1 1 0_1 dc V; down
   R1 1 2 R1; right
   R2 2 0 R2; down
   P1 2_2 0_2; down
   W1 2 2_2; right
   W2 0_1 0; right
   W3 0 0_2; right

Here, P1 defines a port.  This is shown as a pair of open blobs.

Here's the resulting schematic:

.. image:: examples/schematics/voltage-divider.png
   :width: 5cm

Many other components can be drawn than can be simulated.  This
includes non-linear devices such as transistors and diodes and time
varying components such as switches.  For example, here's a common
base amplifier,

.. image:: examples/schematics/common-base.png
   :width: 7cm

This is described by the netlist::

    Q1 3 0 2 pnp; up
    R1 1 2;right
    R2 4 0_4;down
    P1 1 0_1;down
    W 0_1 0;right
    W 0 0_4;right
    W 3 4;right


Network schematics
------------------

One port networks can be drawn with a horizontal layout.  Here's an example:

   >>> from lcapy import R, C, L
   >>> ((R(1) + L(2)) | C(3)).draw()

Here's the result:

.. image:: examples/networks/pickup.png
   :width: 5cm

The s-domain model can be drawn using:

   >>> from lcapy import R, C, L
   >>> ((R(1) + L(2)) | C(3)).smodel().draw()

This produces:

.. image:: examples/networks/pickup-s.png
   :width: 5cm

Internally, Lcapy converts the network to a netlist and then draws the
netlist.  The netlist can be found using the netlist method, for example,

   >>> from lcapy import R, C, L
   >>> print(((R(1) + L(2)) | C(3)).netlist())

yields::

   W 1 3; right, size=0.5
   W 3 4; up, size=0.4
   W 3 5; down, size=0.4
   W 6 2; right, size=0.5
   W 6 7; up, size=0.4
   W 6 8; down, size=0.4
   R 4 9 1; right
   W 9 10; right, size=0.5
   L 10 7 2 0; right
   C 5 8 3 0; right

Note, the components have anonymous identifiers.


IPython Notebooks
=================

IPython notebooks allow interactive markup of python code and text.  A
number of examples are provided in the `lcapy/doc/examples/notebooks`
directory.  Before these notebooks can be viewed in a browser you need
to start an IPython notebook server.

   >>> cd lcapy/doc/examples/notebooks
   >>> ipython notebook

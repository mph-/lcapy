========
Tutorial
========


Introduction
============

mcircuit is a set of circuit analysis classes for Python.  It will
only solve linear, time invariant networks.  In other words, networks
comprised of basic circuit elements (R, L, C, etc.) that do not vary
with time.

It does not support non-linear devices such as diodes or transistors
although it does support simple opamps without saturation.

Mcircuit uses Sympy (symbolic Python) for its values and expressions
and thus the circuit analysis can be performed symbolically.

Internally, the circuit components are stored using their s-domain
equivalents, such as impedances and admittances.  This is convenient
for frequency response analysis but requires an inverse Laplace
transform for transient response analysis.


Preliminaries
=============

- You will need to install mcircuit first.  Currently, this a single
  Python file so you just need to copy it to your working directory or
  set your `PYTHON_PATH` environment variable to find its directory.

- Then fire up your favourite python interpreter, for example, ipython.



Simple circuit element combinations
===================================

Here's an example of resistors in series

   >>> from mcircuit import *
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
method will simplify the netwrok and combine the resistors into a
single resistor equivalent.


Here's an example of a parallel combination of resistors.  Note that
the parallel operator is `|` instead of the usual `||`.

   >>> from mcircuit import *
   >>> Rtot = R(10) | R(5)
   >>> print(Rtot)
   R(10) | R(5)
   >>> print(Rtot.simplify())
   R(10/3)

The result can be performed symbolically, for example,

   >>> from mcircuit import *
   >>> Rtot = R('R_1') | R('R_2')
   >>> print(Rtot)
   R(R_1) | R(R_2)
   >>> print(Rtot.simplify())
   R(R_1*R_2/(R_1 + R_2))

Here's another example using inductors in series

   >>> from mcircuit import *
   >>> L1 = L(10)
   >>> L2 = L(5)
   >>> Ltot = L1 + L2
   >>> print(Ltot)
   L(10) + L(5)
   >>> print(Ltot.simplify())
   L(15)


Finally, here's an example of a parallel combination of capacitors

   >>> from mcircuit import *
   >>> Ctot = C(10) | C(5)
   >>> print(Ctot)
   C(10) | C(5)
   >>> print(Ctot.simplify())
   C(15)


Impedances
==========


Let's consider a series R-L-C network

   >>> from mcircuit import *
   >>> N = R(5) + L(20) + C(10)
   >>> print(N)
   R(5) + L(20) + C(10)
   >>> pprint(N.Z)
        2           
   200⋅s  + 80⋅s + 1
   ─────────────────
          20⋅s    
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


Simple transient analysis
=========================

Let's consider a series R-C network in series with a voltage source

   >>> from mcircuit import *
   >>> N = V(20) + R(5) + C(10)
   >>> print(N)
   V(20) + R(5) + C(10)
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

   >>> from mcircuit import *
   >>> N = V('V_1') + R('R_1') + C('C_1')
   >>> print(N)
   V(V_1) + R(R_1) + C(C_1)
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

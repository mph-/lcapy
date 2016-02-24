=========
Internals
=========


Circuits
========

Circuits are represented using netlists of Mnacpts.  These are
wrappers around Oneports.  Analysis is performed with modified nodal
analysis (MNA).


Networks
========

Networks are comprised of Oneport and Twoport components and are
stored as an abstract syntax tree (AST).

The main attributes are `Voc` (open-circuit s-domain voltage), `Isc`
(short-circuit s-domain voltage), `Z` (s-domain impedance), `Y`
(s-domain admittance).  In addition, there is `Vocac` (open-circuit
phasor voltage) and `Iscac` (short-circuit phasor current).  In
addition, `I` is the component through the one-port (zero by
definition) and `V` is equivalent to the open-circuit voltage `Voc`.

Formerly all Oneport components were either a Thevenin or Norton
component.  As these components were combined (in series or parallel)
a new Thevenin or Norton component was created.  This was efficient
and worked well.  It was also more robust when converting zero
impedances to admittances and vice-versa.  For example, 1 / (1 / 0)
should give 0. However, it is tricky to handle superposition of
multiple independent sources, say ac and dc.  So instead, the same
circuit analysis is performed as for Circuit objects by converting the
network to a netlist.

The current approach needs more thought to deal with the ac/dc/laplace
hackery.  A Component returns s-domain expressions for the `Voc` and
`Isc` methods.  



Values and Expressions
======================

Lcapy uses a number of classes to represent a value or expression.
These classes all inherit from the 'Expr' base class; this is a
wrapper for a sympy expression.  Unfortunately, sympy does not provide
a generic sympy expression class so 'Expr' stores the sympy expression
as its 'expr' attribute.

'cExpr' represents a constant, such as the resistance of a resistor.
The constant must be real and positive.

'tExpr' represents a time domain expression.   This should be real.

'sExpr' represents an s-domain expression.   This can be complex.

'tsExpr' represents an s-domain expression but one that is first
interpreted as a time domain expression.


Note sympy considers that the two expressions x1 and x2 are different
since they have different conditions, where

  x1 = sym.symbols('x')
  x2 = sym.symbols('x', real=True)

Thus x1 - x2 does not simplify to zero.  To overcome this problem,
Lcapy maintains a symbol cache and tries to replace symbols with their
first definition.  The downside is that this may prevent
simplification if the symbol is first defined without any conditions.


C1 1 2 C v0

Here C and v0 should both be real.

V1 1 2 dc V

Here V should be real since there is no phase at DC.

V1 1 2 V

Here V should be real since it should be interpreted as a time domain
voltage.

V1 1 2 s V

Here V can be complex since it is an s-domain voltage.

Things get tricky if we allow a complex s-domain symbol....

Lcapy maintains a set of symbols for each circuit plus a set of
additional symbols defined when creating other objects, such as V
or C.  Symbol names are converted into a canonical format, V1 -> V_1.

sympy defines symbols Q, C, O, S, I, N, E.  It also consider E1 as the
generalized exponential integral function. Thus sympify(5 * E1) fails.

Assumptions are useful for sympy to simplify expressions.  For
example, knowing that a symbol is real or real and positive.

Perhaps symbols cannot be complex or negative?  Complex expressions
could still be created, say Z = X + I * Y or Z = A * exp(phi).

We still have symbols that we know are positive, say a resistance, and
this helps with simplification.  However, initial conditions can be
negative.

C1 1 2 C v0

Here C is real and positive but v0 is real and can be negative.


Assumptions
===========

Assumptions are required to simplify expressions and to help with
inverse Laplace transforms.

There are two types of assumptions.
1. Assumptions used by Sympy, such as real, positive, etc.
2. Assumptions used by Lcapy, such as dc, real, causal, etc.


Sympy assumptions
-----------------

To confuse matters, Sympy has two assumptions mechanisms, old and new.
The old method attaches attributes to symbols, for example,

from sympy import Symbol, Q, exp, I, pi
x = Symbol('x', integer=True)
z = exp(2 * pi * I * x)

The simplify function (or method) uses these attributes.

The new method stores facts, these need not just be about symbols, for
example,

from sympy import Symbol, Q, exp, I, pi
from sympy.assumptions.assume import global_assumptions

x = Symbol('x')
global_assumptions.add(Q.integer(x))
z = exp(2 * pi * I * x)
z = z.refine()

The new method has the advantage that we can collect facts about a
symbol, say from different nets in a netlist.  Since they refer to the
same symbol, there is no problem updating these facts.  The big
problem is how to deal with context, say if we are analysing two
circuits at the same time.  The simplest approach is to create a
context for each circuit and to switch the global_assumptions.

A resistor should have a positive resistance, but what about {a - b}.
We could add an assumption that a - b > 0 but we cannot assume that
both a and b are positive.  Unfortunately, this is the status quo but
is uncommon.


Lcapy assumptions
-----------------

Lcapy expressions have associated assumptions, ac, dc, and causal.
These influence how the result of an inverse Laplace transform is
determined for :math:`t < 0`.

These assumptions are currently not propagated during expression
manipulation.  If so, do we check the assumptions during tests for
equality?

Rather than propagating assumptions, Lcapy assigns them to expressions
after circuit analysis.

Lcapy may need a Superposition class that acts as a container for
different expressions, each with their own assumptions.  For example,
a dc signal plus an ac signal.  Currently, when Vphasor (with ac=True)
and Vs (with dc=True) are added we get a Vs but with no assumptions.


Adding new components
=====================

# Define in grammar.py
# Add class in mnacpts.py for simulation
# Add class in schemcpts.py for drawing


Schematic layout
================

The current layout algorithm assumes that all one-port components such
as resistors and diodes are stretchy.  The x and y positions of
component nodes are determined independently using directed acyclic
graphs.

The steps of the algorithm are:
1. Construct a graph where the edges are the components.  Electrical
   nodes with a common x or y position are combined to reduce the
   graph size.
2. Find longest path through graph.  This determines the maximum
   dimension.  Nodes along this longest path are assigned positions
   based on the maximum distance from the start.  Note, there may be
   multiple parallel paths of the same length; it does not matter
   which is chosen.
3. For each component with an unknown position, find the longest path
   in both forward and backward directions to a node with a known
   position.  This path is traversed counting the number of stretchy
   components and summing their sizes.  Using the distance between the
   positions of the known nodes the stretch per stretchy component can
   be calculated and thus the position of the node.  If the component
   has a dangling node the stretch is zero.


Circuit analysis
================

Lcapy performs circuit analysis using modified nodal analysis (MNA) of
the Laplace representations of signals and impedances.



DC analysis
-----------

If all the independent sources are DC, Lcapy uses some sleight of
hand.  First capacitors are replaced with open-circuits and inductors
are replaced with short-circuits.  Laplace analysis is then performed
but since there are no reactive components, there are no transients
produced.  Thus all the voltages and currents have the form :math:`X /
s`.  where :math:`X` is does not depend on :math:`s`.  Finally, when
taking an inverse Laplace transform, the result is extrapolated for
:math:`t < 0`.


AC analysis
-----------

AC analysis uses phasors to represent voltages and currents and
impedances to avoid differential equations.  Using Laplace analysis
for AC analysis is tricky since we need to determine initial values
since Laplace analysis assumes that the sources are zero for :math:`t
< 0`.  Without the initial values, Laplace analysis can generate bogus
transients.

When Lcapy detects an AC circuit analysis problem, it uses phasors of
angular frequency :math:`\omega_1`.  Note, :math:`\omega` is the
angular frequency Lcapy uses for the Fourier domain.  Thus the `.V`
and `.I` attributes return phasors for AC problems and s-domain
signals for Laplace problems.


General analysis
----------------

A circuit with a variety of independent sources, for example, ac and
dc, the result can be determined using superposition of each source
treated independently.  Currently, Lcapy does not do this
automatically.  Instead, it is necessary to use the `kill_except`
method to kill the other sources.

General analysis is tricky with arbitrary switching sources since each
change of amplitude, phase, or frequency generates a transient.
Analysis is straightforward if all sources only switch once but if
they switch multiple times we need to consider a number of dependent
initial value problems.



Expression manipulation
-----------------------

cos(x).rewrite(exp) ->  exp(j*x) / 2 + exp(-j*x)/2
(exp(j*x) / 2 + exp(-j*x)/2).rewrite(cos) -> cos(x)
(exp(j*x) / 2 + exp(-j*x)/2).rewrite(sin) -> cos(x)

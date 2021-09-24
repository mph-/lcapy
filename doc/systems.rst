.. _systems:

=======
Systems
=======


State-space representation
==========================

A state-space object is created from the state matrix, `A`, input
matrix, `B`, output matrix `C`, and feed-through matrix `D`::

    >>> ss = StateSpace(A, B, C, D)

A state-space object can also be created from lists of the numerator
and denominator coefficients `b` and `a`::

   >>> ss = StateSpace.from_ba(b, a)

By default, the controllable canonical form CCF is created.  The
observable canonical form OCF is created with::
  
   >>> ss = StateSpace.from_ba(b, a, form='OCF')
   

State-space from transfer function
----------------------------------

Transfer functions (and impedances and admittances) can be converted
to a state-space representation.  Here's an example::

   >>> Z = (s**2 + a) / (s**3 + b * s + c)
   >>> ss = Z.state_space('CCF')
   
State-space representation are not unique; Lcapy uses the controllable
canonical form (CCF) and the observable canonical form (OCF).  The CCF
form of the state-space matrices are::

   >>> ss.A
   ⎡0   1   0⎤
   ⎢         ⎥
   ⎢0   0   1⎥
   ⎢         ⎥
   ⎣-c  -b  0⎦

   >>> ss.B
   ⎡0⎤
   ⎢ ⎥
   ⎢0⎥
   ⎢ ⎥
   ⎣1⎦

   >>> ss.C
   [a  0  1]

   >>> ss.D
   [0]


Transfer function from state-space
----------------------------------

For a single-input single-output (SISO) system the transfer function
is obtained with the `transfer_function()` method, for example::

   >>> ss = StateSpace(A, B, C, D)
   >>> G = ss.transfer_function
   

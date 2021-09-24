.. _systems:

=======
Systems
=======


State-space representation
==========================

Transfer functions (and impedances and admittances) can be converted
to a state-space representation.  Here's an example::

   >>> Z = (s**2 + a) / (s**3 + b * s + c)
   >>> ss = Z.state_space('CCF')
   
State-space representation are not unique; this uses the canonical
controllable form (CCF).  With this form the state-space matrices
are::

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

   
   




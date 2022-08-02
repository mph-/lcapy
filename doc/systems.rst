.. _systems:


=======
Systems
=======

.. _state-space:


State-space
===========


Lcapy has two state-space representations: `StateSpace` for
continuous-time linear time-invariant systems and `DTStateSpace` for
discrete-time linear time-invariant systems.  Both representations
share many methods and attributes.


Continuous-time state-space representation
==========================================

A state-space object is created from the state matrix, `A`, input
matrix, `B`, output matrix `C`, and feed-through matrix `D`::

    >>> ss = StateSpace(A, B, C, D)

A state-space object can also be created from lists of the numerator
and denominator coefficients `b` and `a`::

   >>> ss = StateSpace.from_ba(b, a)

By default, the controllable canonical form CCF is created.  The
observable canonical form OCF is created with::

   >>> ss = StateSpace.from_ba(b, a, form='OCF')

Similarly, the diagonal canonical form DCF is created with::

   >>> ss = StateSpace.from_ba(b, a, form='DCF')

For the DCF, the poles of the transfer function must be unique.


State-space from transfer function
----------------------------------

Transfer functions (and impedances and admittances) can be converted
to a state-space representation.  Here's an example::

   >>> Z = (s**2 + a) / (s**3 + b * s + c)
   >>> ss = Z.state_space('CCF')

State-space representation are not unique; Lcapy uses the controllable
canonical form (CCF), the observable canonical form (OCF), and the
diagonal canonical form (DCF).  The CCF form of the state-space
matrices are::

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


Discrete-time state-space representation
========================================

Discrete-time state-space objects are defined in a similar manner to
continuous-time state-space objects and share many methods and
attributes.  A discrete-time state-space object is created from the
state matrix, `A`, input matrix, `B`, output matrix `C`, and
feed-through matrix `D`::

    >>> ss = DTStateSpace(A, B, C, D)

A state-space object can also be created from lists of the numerator
and denominator coefficients `b` and `a`::

   >>> ss = DTStateSpace.from_ba(b, a)

By default, the controllable canonical form CCF is created.  The
observable canonical form OCF is created with::

   >>> ss = DTStateSpace.from_ba(b, a, form='OCF')

Similarly, the diagonal canonical form DCF is created with::

   >>> ss = DTStateSpace.from_ba(b, a, form='DCF')

For the DCF, the poles of the transfer function must be unique.


For example::

   >>> ss = DTStateSpace(((0, 1), (1, 0)), (1, 1), ((1, 2), ), [1])

   >>> ss.A
   ⎡0  1⎤
   ⎢    ⎥
   ⎣1  0⎦

   >>> ss.B
   ⎡1⎤
   ⎢ ⎥
   ⎣1⎦

   >>> ss.C
   [1  2]

   >>> ss.D
   [1]

   >>> ss.state_equations()
   ⎡x₀(n + 1)⎤   ⎡0  1⎤ ⎡x₀(n)⎤   ⎡1⎤
   ⎢         ⎥ = ⎢    ⎥⋅⎢     ⎥ + ⎢ ⎥⋅[u₀(n)]
   ⎣x₁(n + 1)⎦   ⎣1  0⎦ ⎣x₁(n)⎦   ⎣1⎦

   >>> ss.output_equations()
                    ⎡x₀(n)⎤
   [y₀(n)] = [1  2]⋅⎢     ⎥ + [1]⋅[u₀(n)]
                    ⎣x₁(n)⎦



   >>> ss.controllability_matrix
   ⎡1  1⎤
   ⎢    ⎥
   ⎣1  1⎦

   >>> ss.is_controllable
   False


   >>> ss = DTStateSpace(((0, 1), (1, 1)), (1, 1), ((1, 2), ), [1])

   >>> ss.A
   ⎡0  1⎤
   ⎢    ⎥
   ⎣1  1⎦

   >>> ss.B
   ⎡1⎤
   ⎢ ⎥
   ⎣1⎦

   >>> ss.C
   [1  2]

   >>> ss.D
   [1]

   >>> ss.is_stable
   False

   >>> ss.eigenvalues
   [-1, 1]

   >>> ss.controllability_matrix
   ⎡1  1⎤
   ⎢    ⎥
   ⎣1  2⎦

   >>> ss.is_controllable
   True

   >>> ss.is_observable
   True

   >>> ss.state_transfer([[2], [3]], xinitial=[0, 0])
   ⎡5⎤
   ⎢ ⎥
   ⎣7⎦

   >>> ss.minimum_energy_input(2, [5, 7], [0, 0])
   ⎡2⎤
   ⎢ ⎥
   ⎣3⎦

   >>> ss.minimum_energy(2, [5, 7], [0, 0])

   >>> ss.minimum_energy_input(3, [5, 7], [0, 0])
   ⎡5/3⎤
   ⎢   ⎥
   ⎢1/3⎥
   ⎢   ⎥
   ⎣4/3⎦

   >>> ss.minimum_energy(3, [5, 7], [0, 0])
   14/3


State-space operations
======================


Model balancing
---------------

This returns a new StateSpace object that has the controllability and
observability gramians equal to the diagonal matrix with the
Hankel singular values on the diagonal.  For example::

   >>> ss2 = ss.balance()

Note, this requires numerical A, B, C, D matrices.


Model reduction
---------------

A balanced reduction can be performed using::

   >>> ss2 = ss.balance_reduce(threshold=0.1)

where states are removed with a Hankel singular value below the
threshold.   Note, this requires numerical A, B, C, D matrices.

Alternatively, specific states can be removed.  For example::

     >>> ss2 = ss.reduce(elim_states=[1, 3])

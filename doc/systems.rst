.. _systems:


=======================
Continuous-time systems
=======================


Linear time invariant filters
=============================

Linear time invariant filters are created with the `LTIFilter` class::

    >>> fil = LTIFilter(b, a)

where `b` is a list or array of the transfer function numerator
 coefficients and `a` is a list or array of the transfer function
 denominator coefficients.  For example,

    >>> fil = LTIFilter(('b0', ), ('a0', 'a1'))

The coefficients can be found from the `a` and `b` attributes::

   >>> fil.a
   (a₀, a₁)

   >>> fil.b
   (b₀,)

The filter's transfer function is found with the `transfer_function()` method::

   >>> fil.transfer_function()

       b₀
   ─────────
   a₀⋅s + a₁

its frequency response is found with the `frequency_response()` method::

   >>> fil.frequency_response()
          b₀
   ───────────────
   2⋅ⅉ⋅π⋅a₀⋅f + a₁

its impulse response is found with the `impulse_response()` method::

   >>> fil.impulse_response()
       -a₁⋅t
       ──────
         a₀
   b₀⋅ℯ      ⋅u(t)
   ───────────────
          a₀

its step response is found with the `step_response()` method::

   >>> fil.step_response()
      ⎛      -a₁⋅t ⎞
      ⎜      ──────⎟
      ⎜        a₀  ⎟
      ⎜1    ℯ      ⎟
   b₀⋅⎜── - ───────⎟⋅u(t)
      ⎝a₁      a₁  ⎠

The filter's differential equation is found with the `differential_equation()` method::

   >>> fil.differential_equation()
                  d
   a₀⋅y(t) = - a₁⋅──(y(t)) + b₀⋅x(t)
                  dt
The input and output symbols can be changed with the `inputsym` and `outputsym` arguments.

The response due to intial conditions is found with the `initial_response()` method:

   >>> fil.initial_response(('y0', ))
          ⎛           -a₁⋅t      ⎞
          ⎜           ──────     ⎟
          ⎜             a₀       ⎟
          ⎜δ(t)   a₁⋅ℯ      ⋅u(t)⎟
   -a₁⋅y₀⋅⎜──── - ───────────────⎟
          ⎜ a₀            2      ⎟
          ⎝             a₀       ⎠


Continuous-time linear time invariant filter attributes
-------------------------------------------------------

- `a` denominator coefficients as a list
- `b` numerator coefficients as a list
- `is_marginally_stable` True if impulse response marginally stable
- `is_stable` True if impulse response stable


Continuous-time linear time invariant filter methods
----------------------------------------------------

- `differential_equation()` creates discrete-time differential equation
- `impulse_response()` creates discrete-time domain impulse response
- `initial_response()` returns discrete time-domain response due to initial conditions
- `response()` returns discrete time-domain response due to input signal and initial conditions
- `transfer_function()` creates s-domain transfer function
- `sdomain_initial_response()` returns s-domain response due to initial conditions



Differential equations
======================

Differential equations are represented by the `DifferentialEquation`
class.  They are usually created by the `differential_equation()`
method of the `LTIFilter` class or from transfer functions.  For example::

   >>> fil = LTIFilter(('b0', ), ('a0', 'a1'))
   >>> de = fil.differential_equation()
   >>> de
                  d
   a₀⋅y(t) = - a₁⋅──(y(t)) + b₀⋅x(t)
                  dt

There are two attributes: `lhs` for the left-hand-side and `rhs` for
the right-hand-side,

   >>> de.lhs
   a₀⋅y(t)
   >>> de.rhs
        d
   - a₁⋅──(y(t)) + b₀⋅x(t)
        dt

A transfer function is created with the `transfer_function()` method::

  >>> de.transfer_function()
      b₀
   ─────────
   a₀ + a₁⋅s

An `LTFilter` object is created with the `lti_filter()` method::

  >>> fil = de.lti_filter()


Differential equation attributes
--------------------------------

- `lhs` left-hand-side of the equation
- `rhs` right-hand-side of the equation
- `inputsym` input symbol, usually 'x'
- `outputsym` input symbol, usually 'y'


Differential equation methods
-----------------------------

- `dlti_filter()` creates continuous-time linear time invariant filter (`LTIFilter`) object
- `separate()` separates the input expressions from the output expressions.
- `transfer_function()` creates s-domain transfer function


.. _state-space:


Continuous-time state-space representation
==========================================

Lcapy has two state-space representations: `StateSpace` for
continuous-time linear time-invariant systems and `DTStateSpace` for
discrete-time linear time-invariant systems.  Both representations
share many methods and attributes.

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


=====================
Discrete-time systems
=====================


Difference equations
====================

Difference equations can be generated from transfer functions and
impulse responses.  Both FIR and IIR (direct form I) can be generated.
For example::

  >>> H = (z + 2) / z**2
  >>> H.difference_equation('x', 'y', 'fir')
  y(n) = 2⋅x(n - 2) + x(n - 1)

Difference equations can be created explicitly, for example::

  >>> de = difference_equation('y(n)', '2 * x(n - 2) + x(n - 1)')

The `separate()` method separates the input expressions from the
output expressions.   For example::

  >>> de = difference_equation('y(n)', '2 * y(n - 1) + x(n)')
  >>> de.separate()
  y(n) - 2⋅y(n - 1) = x(n)


Difference equation attributes
------------------------------

- `lhs` left-hand-side of the equation
- `rhs` right-hand-side of the equation
- `inputsym` input symbol, usually 'x'
- `outputsym` input symbol, usually 'y'


Difference equation methods
---------------------------

- `dlti_filter()` creates discrete-time linear time invariant filter (`DLTIFilter`) object
- `separate()` separates the input expressions from the output expressions.
- `transfer_function()` creates z-domain transfer function


Discrete-time transfer functions
================================

A discrete-time transfer functions can be determined from a difference
equation or a DLTI filter.  For example::

   >>> de = difference_equation('y(n)', '2 * x(n - 2) + x(n - 1)')
   >>> H = de.transfer_function()
   >>> H
   z + 2
   ─────
     2
    z


Discrete-time transfer function methods
---------------------------------------

- `dlti_filter()` creates discrete-time linear time invariant filter (`DLTIFilter`) object
- `difference_equation()` creates discrete-time difference equation


.. _DLTIfilter:

Discrete-time linear time invariant filters
===========================================

A discrete-time linear time invariant filter can be specified by its
numerator and denominator coefficients.  For example, a first-order,
discrete-time, recursive low-pass filter can be created with:

   >>> a = symbol('a')
   >>> lpf = DLTIFilter((1 - a, ), (1, -a))

The difference equation can be printed using::

   >>> lpf.difference_equation()
   y(n) = a⋅y(n - 1) + (1 - a)⋅x(n)

The transfer function can be printed using::

   >>> lpf.transfer_function()
   z⋅(a - 1)
   ─────────
     a - z

The impulse response can be printed using::

   >>> lpf.impulse_response()
    n
   a ⋅(1 - a)⋅u[n]

The general response to an input `x(n)` can be printed using::

  >>> lpf.response(x, ni=(0, 5))

For a recursive filter, the initial conditions can also be specified::

  >>> lpf.response(x, ic=[1], ni=(0, 5))

The input to the filter can be a `DiscreteTimeDomainExpression` or a sequence.
The output is a sequence.

A discrete-time LTI filter can be created from difference equations
and transfer functions.   For example::

  >>> de = DifferenceEquation('2 * y(n)', '4 * y(n + 1) - 3 * y(n-3) -2 * x(n) - 5 * x(n-3)')
  >>> fil = de.dlti_filter()
  >>> fil.a
  [4, -2, 0, 0, -3]
  >>> fil.b
  [0, 2, 0, 0, 5]
  >>> fil.difference_equation()


Discrete-time linear time invariant filter attributes
-----------------------------------------------------

- `a` denominator coefficients as a list
- `b` numerator coefficients as a list
- `is_marginally_stable` True if impulse response marginally stable
- `is_stable` True if impulse response stable


Discrete-time linear time invariant filter methods
--------------------------------------------------

- `difference_equation()` creates discrete-time difference equation
- `impulse_response()` creates discrete-time domain impulse response
- `initial_response()` returns discrete time-domain response due to initial conditions
- `inverse()` creates an inverse filter by switching numerator and denominator coefficients
- `response()` returns discrete time-domain response due to input signal and initial conditions
- `transfer_function()` creates z-domain transfer function
- `zdomain_initial_response()` returns z-domain response due to initial conditions



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

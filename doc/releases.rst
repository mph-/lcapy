=============
Release Notes
=============

V 0.79
======

- Fix units for `delta(x)`, `diff(f, x)`, `integrate(f, x)`.

- `state.canonical_form` controls whether units are stored in canonical form, e.g., watt rather than volt * ampere.



V 0.78
======

Tracking and checking of units for quantities is functional, for example::

   >>> state.show_units = True
   >>> V = voltage(4)
   >>> Z = impedance(2)
   >>> I = V / Z
   >>> I
   2⋅A
   >>> state.abbreviate_units = False
   >>> I
   2⋅ampere
   >>> I.units
   ampere

You cannot add two expressions with different units, `current(2) + voltage(4)` will fail.  If `loose_units` is defined (default), then constants can be added to expressions, for example::

  >>> voltage(4) + 1
  5⋅V
  >>> state.loose_units = False
  >>> voltage(4) + 1
  ValueError: Cannot determine ConstantTimeDomainVoltage(4*V) + ConstantDomainExpression(1) since the units V are incompatible with 1

Units are not correctly tracked for function calls, for example, `log(voltage(10)` or `delta(t)`.
   

=============
Release notes
=============

V0.79
=====

- Fix units for `delta(x)`, `diff(f, x)`, `integrate(f, x)`.

- `state.canonical_form` controls whether units are printed in canonical form, e.g., watt rather than volt * ampere.

- `dc`, `ac`, `causal` attributes removed, instead use `is_dc`, `is_ac`, and `is_causal`.

- `dc` returns dc component, `ac` returns ac components as dictionary (this may change), `transient` returns transient component

- Fix expression printing with units if have no units

- Fix expression printing with units if expression is 1

- Improved Laplace transforms for convolutions  


V0.78
=====

- Tracking, checking, and printing units for quantities is functional, for example::

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

- Prevent addition/subtraction of two expressions with different units, `current(2) + voltage(4)` will fail.  If `loose_units` is defined (default), then constants can be added to expressions, for example::

   >>> voltage(4) + 1
   5⋅V
   >>> state.loose_units = False
   >>> voltage(4) + 1
   ValueError: Cannot determine ConstantTimeDomainVoltage(4*V) + ConstantDomainExpression(1) since the units V are incompatible with 1

Units are not correctly tracked for function calls, for example, `log(voltage(10)` or `delta(t)`.
   

Older versions
==============

- V0.77 reverts phase as a quantity and fixes plots.  Component attributes are renamed for consistency (is_resistor etc.).  omega0 is now positive.  Allow Z / Z and Y / Y.  Fix matrices.  Lazily create expression classes.  More unit tests! 

- V0.76 fixes the units and adds many more tests.  Adds phase quantity.  Fixes phasors.

- V0.75 introduces a major change to expression classes with tighter restrictions on operations between expressions.  For example, a current expression cannot be added to a voltage expression.  There is also experimental support for showing units.  Added phasor domain.  Discrete-time support is now enabled.  This introduces three new domain variables, n, k, and z.  More Fourier transforms added.  Sinc and rect functions added.

- V0.74 supports simplification of netlists, adds more rigorous type checking for expressions, improve printing of conditional expressions.

- V0.73 improves printing of Voltage and Current, adds phasor attributes to Voltage and Current, fixes magnitude and phase for Phasor, fixes printing of greek symbols, tidies canonical representation, wraps R, X, B, G attributes for Immittance, doc reorganisation.

- V0.72 uses CI for docs plus many assorted bug fixes.

- V0.71 uses much faster matrix inversion (if sympy-1.8 installed) otherwise falls back on ADJ method  instead of the GE method which has a serious time regression with sympy-1.6.2

- V0.70 adds improved nodal and mesh analysis.

- V0.69 adds common-mode gain for opamps and polyphase-twoports.

- V0.67 adds time-stepping simulation, supernode detection, and polyphase circuits.

- V0.66 tidies up two-port parameters.  S and T parameters are
  added.  The A, B, G, H, Y, Z parameters are renamed to Aparams, etc. to avoid conflict with
  matrix transpose and Hermitian transpose operators.  issymmetrical, isshunt renamed to is_symmetrical,
  is+shunt, etc.   Eq, MatMul, MatAdd, Mul, Add functions added.  Expr.__getattr_ converts lists to ExprList.
  Adds symbols attribute to Matrix.  Ensures symbols in immitance default to complex.

- V0.65 introduces random networks.  Adds simplification for DiracDelta and Heaviside.  Adds node checking for Netlist methods.

- V0.64 adds wye-delta, delta-wye transformations.  Adds resistive companion models.  Fix state-space if have no sources.  Fixes assumption merging.  Adds verbatim argument for laplace_transform.   Simplifies mutual inductance.

- V0.63 fixes mirroring of opamps in schematics and introduces mirrorinputs option

- V0.62 adds search, save, annotate_voltage, annotate_current, kill_zero methods.  Fixes solve.

- V0.61 improves Laplace and z-transforms.

- V0.60 replaces DiracDelta with UnitImpulse and Heaviside with UnitStep for discrete-time expressions.

- V0.52 improves the component positioning algorithm for schematics.

- V0.51 improves the domain transformation infrastructure,

- V0.50 changes phasors to have a default angular frequeny of omega_0 instead of omega to avoid confusion with angular frequency in Fourier transforms, adds preliminary phasor plots, improves noise signal classes, improves the infrastructure, and fixes many bugs.

- V0.49 adds mechanical components, better parameterization, faster partical fraction expansion, improved Z transforms, IIR difference equations, and differential equations.

- V0.48 fixes z-transforms, adds better caching for Laplace and z-transforms, convert rational numbers to floats on schematics, fixes expr rpow.

- V0.47 introduces subs method for netlists, initialize method of netlists, better clarification for external programs, removes Y and Z methods for Circuits, removes anon ids from circuit components, adds remove_condition, force_causal, is_conditional, is_rational_function, is_strictly_proper, adds isoamp, inamp, and bug fixes

- V0.42 bug fixes for discrete-time signals

- V0.41 introduces experimental discrete-time signals

- V0.40 fixes schematics

- V0.39 miscellaneous bug fixes

- V0.38 reverts the experimental behaviour of 0.37.  Instead it introduces new classes for general immitances that tries to display them in the most suitable format.

- V0.37 changes the API for admittances and impedances.  The
  attributes Y and Z return the impedance in terms of omega rather
  than s as in the previous versions.  The old behaviour is provided
  with the Ys and Zs attributes (generalized admittance and
  impedance).  It also has better distinction between the impedance of
  a component and the driving point impedance.

- V0.36 has improved handling of complex conjugate poles.        

- V0.34 switched to using setuptools and pushed to https::pypi.org

- V0.33 reworked expression printing infrastructure

- V0.32.3 introduces state-space analysis.  The API is experimental and may change.

- V0.32.0 changes the naming of symbolic values.  Previously R1 was converted to R_1 before being converted into a SymPy symbol.  This behaviour was not obvious for symbol substitution.  Now the symbol names are converted on printing.

- V0.31.0 reworks schematic drawing.  The syntax for chips has changed since there are no explicit nodes in the netlist.

- V0.30.0 tweaks the syntax to perform transformations based on the argument, e.g., V(s) or V(t)

- V0.28.0 works with Sympy 1.2.

- V0.26.0 adds noise analysis.

- V0.25.1 adds time-domain analysis for circuits without reactive components.

- From version 0.25.0, Lcapy performs more comprehensive circuit analysis using combinations of DC, AC, and Laplace analysis.  This added functionality has resulted in a slight change of syntax.  cct.R1.V no longer prints the s-domain expression but the decomposition of a signal into each of the transform domains.


Future plans
============

- As of V0.79 there are no envisaged modifications to the API :-)

- Add more unit tests

- Add more Laplace and Fourier transforms
  


  

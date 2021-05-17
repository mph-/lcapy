=============
Release notes
=============

V0.90
=====

- Adds call notation to access element of `Sequence`

- Adds `as_array()` method for `Sequence`

- Modifies `evaluate()` method for `Sequence`  

- Adds DTFilter

- Adds evalf method for container classes
 
- Fixes access of element in a sequence

- Adds override argument to expr

- User defined symbols override SymPy symbols
  
- Does not print user defined symbols in canonical form

- Reworks equation function

- Removes undefs when simplifying or solving

- Fixes inverse z-transforms for z**n

- Adds many new z-transforms

V0.89
=====

- Adds title arg for plots

- Fixes label args for pole zero plots

- Adds periodic sinc function
  
- Adds normalised and unnormalised versions of sinc

- Fixes evaluation of sinc

- Fixes phasors with no var  
  

V0.88
=====

- Evaluates unit step

- Adds new z-transforms

- Fixes inverse z-transform of repeated pole

- Ensures discrete-time string conversions converted

- Adds `tri(t)` and `trap(t, alpha)` functions

- Adds new Fourier transforms  

- Fixes `(rect(t) * cos(2 * pi * t))(f)`

- Fixes `rect(t)(f)`

- Functions return `Expr` objects

  
V0.87
=====

- Fixes general problems with phasor transforms

- Adds `bode_plot()` method to phasors and s-domain expressions

- Adds `pole_zero_plot()` method to s-domain expressions  

- Allows complex signals to be considered as ac signals  

- Adds `is_complex_signal` attribute

- Documents transformations

- Allows `sexpr(voltage(4))` as well as `voltage(sexpr(4))`, etc.

- Add `links()` method to `CircuitGraph`


V0.86
=====

- Enables short-cut for transforming s to jw or w domains

- Adds noiseless resistors

- Adds subs() method for networks

- Adds noisy() method for networks  

- Adds T arg to noisy() methods
  

V0.85
=====

- Supports SymPy 1.8

- Changes behaviour of V1 1 2 to be equivalent to V1 1 2 V1.  The same
  applies for I1 1 2.  This is consistent with other component
  definitions and allows netlist substitutions.

- Allows substitutions for constant expressions  

- Fixes is_unchanging for phasors

- Adds additional opamp noise tutorials

- Fixes frequency plots

- Reworks `CircuitGraph` to suport trees

- Changes `CircuitGraph` `nodes()` method to be an attribute 

- Fixes state-space analysis with current source

- Adds differential drivers to schematics

- Adds `has()` and `replace()` methods to netlists

- Allows component names to specified as well as nodes for the `transfer()` method


V0.84
=====

- Adds debugging support when generating schematics
  
- Reverts to using temporary dictionary for temporary files during schematic generation

- Ensures log file closed before deleting

- Fixes units for 1/s  


V0.83
=====

- Adds new opamp tutorials on transimpedance amplifiers and multi-feedback filters

- Adds an experimental component placement algorithm for schematics

- Schematics are converted to pdf in the local directory to access relative files

- Adds support for PGF files to be included into schematics with the image keyword

- Improves some Laplace transforms  

- Fixes state-space model for current sources


V0.82
=====

This release primarily improves the component placement algorithm for schematics that also prevents crashes

- Improve component placement algorithm; add message suggesting constraint component to ensure symmetry

- Improve component placement graphs for debugging

- Require pdflatex for schematic tests
  

V0.81
=====

This is mostly bug fixes

- Add tests for loop and nodal analysis

- Add tests for schematics  

- Improve twoport printing


V0.80
=====

This is mostly bug fixes

- Require sympy > 1.7.1

- Install ghostscript for tests  

- Fix IDFT X(k)

- Add tests for CircuitGraph

- Simplify products of u(t)  

- Add tests for sinc, rect 

- Fix convolution units

- Fix FT of convolution


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
  


  

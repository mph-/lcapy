=============
Release notes
=============

V1.12
=====

- Adds support for lcapy-tk (this is a GUI under development for drawing and analyzing schematics with Lcapy, see https://github.com/mph-/lcapy-gui)

- Adds connection attributes to annotations

- Adds attribute definitions (see :ref:`attribute_definitions`)


V1.11
=====

- Simplifies Circuitikz output for schematics

- Allows autonaming for netlist components (see :ref:`autonaming`)

- Schematic attributes specified in the last netlist entry are considered first

- Use lower case v for time-domain nodal analysis

- Use lower case i for time-domain loop analysis

- Fixes loop analysis

- Adds new schematic syntax for drawing nodes and implicit connections (see :ref:`node_attributes`)


V1.10.1
=======

 - Works with NumPy 1.24.0

 - Fixes MOSFET drawing


V1.10
=====

- Adds reluctance component RL for drawing

- Adds parameter estimation method `estimate()` to expressions (see
  :ref:`parameter_estimation`)

- Disables png output for Jupyter

- Unify `resistance()`, `conductance()`, `capacitance()`,
  `inductance()`, `susceptance()`, and `reactance()` to return
  `ConstantFrequencyResponseDomain` objects

- Fixes units after integration and differentiation

- Updates printing to be compatible with SymPy printing API changes

- Modifies scaling for discretization of continuous-time signals (for admittance, impedance, and transfer function quantities)

- Fixes stability checks

- Adds `dlti_filter()` method for time-domain expressions


V1.9
====

- Supports other iterables for `subs()`

- Adds `approximate_dominant()` method to expressions (see
  :ref:`approximation`)

- Adds units to parameterization definitions

- Adds units to poles and zeros

- Adds blocks for schematics

- Adds tunable kind for schematic components

- Adds chokes for schematics

- Introduces phasor ratio, frequency response, and angular frequency response domains

- Changes `jw` to be the domain variable for the angular frequency response domain

- Adds `jf` domain variable for the frequency response domain

- Warns if old version of Circuitikz found


V1.8
====

- Compatible with SymPy-1.11

- Converts s * t, f * t, w * t to time domain with warning

- Makes phasor arithmetic stricter

- Fixes phasor ratios

- Fixes Bode plot of phasors

- Adds var argument to `bode_plot()` for linear/angular frequency

- Adds Nichols plot

- Adds `j2pif`

- Use `frequency_response()` method for Bode plots (this does not generate Dirac deltas for marginally stable systems)

- Simplifies magnitude of expression with Dirac delta terms

- Fixes plotting of expressions with Dirac deltas outside desired region

- Fixes Laplace to Fourier shortcut

- Adds `is_marginally_stable` attribute to expressions

- Adds `remove_disconnected()`, `remove_dangling()`, `remove_dangling_wires()`

- Adds `select` and `ignore` argument to `simplify()`, `simplify_series()`, `simplify_parallel()`

- Adds `is_dangling` and `is_disconnected` attributes to components

- Warns if using `I` for current source value (this is considered the imaginary operator by SymPy)

- Ensures unique names chosen

- Reduces recursion depth when trying to draw bogus schematics

- Removes checks for ubuntu-18.04


V1.7
====

- Adds `convert_IVP()` method to convert a circuit with switches to an initial value problem

- Handles DC analysis for capacitors by adding a conductance in parallel and considering the limit as the conductance goes to zero

- Adds `replace_switches()` and `replace_switches_before()` to remove switches from a circuit

- Adds `switching_times()` to determine the times when switches activate

- Fixes `mirror` and `invert` attributes for SPDT switches

- Improves debugging for conversion of schematic to png

- Fixes lower limit of convolution when using ILT

- Adds comparison for equations

- Fixes z-domain frequency response

- Adds `LTIFilter` and `DifferentialEquation` classes

- Fixes definition of `psinc()`

- Adds `abc` module to mimim SymPy

- Inherits functions docs from SymPy


V1.6
====

- Fixes autoground for nodes that are not drawn (e.g., with opamp)

- Fixes solving system of equations in Laplace domain

- No longer assumes zero initial conditions for Laplace transforms of
  derivatives

- Adds `zero_initial_conditions` argument for Laplace transforms

- Adds `limit` function

- Fixes initial conditions for loop and nodal analysis

- Fixes `U`, `X`, and `X0` attributes for state space analysis


V1.5.1
======

- Fixes drawing of implicit nodes

- Adds node_label_anchor for repositioning of node labels


V1.5
====

- Uses SymPy-1.10.1 with improved Laplace transform support

- Adds implicit connections for oneport components in netlists, see :ref:`implicit_connections`

- Adds autoground for schematics, see :ref:`autoground`

- Improves choice of node names for nodal analysis

- Avoids double subscripts for LaTeX output

- Adds named parameters for netlists, such as `E1 1 0 opamp 2 3 Ro=Ro`

- Models fully differential and instrumentation amplifiers

- Modifies transistor sizes and improve transistor labelling to work around Circuitikz changes

- Improves math-mode detection for labels

- Adds `0V` implicit connection

- Tidies naming on schematics if the value is the same as the name

- Adds `degrees` and `radians` functions

- Adds `nsolve()` method for numerical solving

- Increases dpi for schematics to 300

- Adds more Fourier transforms for functions of exponentials

- Adds `is_stable` and `is_realizable` attributes

- Unwraps phase for Bode plots

- Removes `omega0` from domain variables

- Ignores `ac` and `dc` assumptions for inverse Laplace transforms

- Adds `kill_noise()` method

- Ignores small imaginary part for `fval` and warns about larger imaginary parts

- Fixes phasor decomposition

- Ensures real symbols are positive by default

- Adds `kind` attribute to voltage/current sources


V1.4
====

- `color` attribute applies to whole schematic; use `help_lines_color` to specify the color of the help lines

- `in_series` and `in_parallel` return lists rather than sets

- Fixes node renumbering when have chips

- Adds `annotate()` method for circuits

- Warns about matrix inversion time for large matrices

- Warns about degenerate circuits

- Fixes state-space analysis when there are no state variables

- Renames `short` to `short-circuit` and adds `open-circuit`

- Adds `voltage_gain()`, `current_gain()`, `transadmittance()`, `transimpedance()` methods for netlists

- Adds `voltage_gain`, `current_gain`, `transadmittance`,
  `transimpedance`, `forward_forward_voltage_gain`,
  `forward_current_gain`, `forward_transadmittance`,
  `forward_transimpedance`, `reverse_voltage_gain`,
  `reverse_current_gain`, `reverse_transadmittance`,
  `reverse_transimpedance`  attributes for networks

- Adds `apply_test_current()` and `apply_test_voltage()` methods

- Fixes `voltage_dir` argument for schematics

- Adds symbol registry

- Shares symbol registry for all circuits

- Allows fancy symbol names

- Checks if components connected if MNA fails

- Adds `wired_to` and `is_wired_to` attributes

- Fixes `nosim` argument for diodes and transistors

- Adds `TLlossless` for lossless transmission lines

- Adds transient response at start of transmission line


V1.3
====

- Adds support for more transistor types in schematics

- Warns if there are no sources in circuit analysis

- Warns if use `k` for coupling coefficient

- Fixes force option for `symbol()`

- Adds Laplace transforms for `ramp`, `rampstep`, `rect`, `tri`

- Adds `ramp()` and `rampstep()` functions

- Adds `expand_functions()` method to `Expr`

- Renames `expandresponse()` to `expand_response`

- Fixes setting causal assumption when extracting from a superposition

- Adds `plot_deltas` argument to `plot()` methods

- Avoids wrapping Jupyter notebook result

- Adds preliminary support for triodes

- Tries harder to find poles and zeros

- Improves finding numerator and denominator expressions

- Fixes conversion to norm Fourier and norm angular Fourier domains

- Makes result of difference equation causal

- Fixes `transfer_function()` and `impulse_response()` for `DLTIFilter`

- Fixes Z-transform for down-sampling

- Fixes discrete-time convolution

- Allows `(f)` notation for DTFT

- Adds lossless transmission line component

- Adds `short()` method to `Circuit`

- Adds `in_series()` and `in_parallel()` methods for components


V1.2.4
======

- Lazily import scipy, numpy, and networkx to speed up loading

- Allows two-ports to be created from netlist using component names


V1.2.3
======

- Fixes voltage and current source drawing for new CircuiTikz

- Adds inverse Laplace transforms for lossless transmission line responses

- Adds `nosim` attribute to ignore component in analysis

- Warns if current name is I


V1.2.2
======

- Adds inverse Laplace transforms for reciprocals of hyperbolic functions

- Fixes printing of reasons for MNA failure

- Fixes `ignore` attribute for schematics

- Renames `TxLine` to `TransmissionLine`

- Adds Z-transform for down-sampling

- Applies similarity and shift theorems for Fourier transforms

- Determines roots numerically if cannot be found symbolically

- Fixes default plot type for frequency plots

- Adds `MatMul` and `MatAdd` functions

- Adds `Z1sc`, `Z2sc`, `Z1oc`, `Z2sc`, etc. for each two-port model

- Adds `Transformer` two-port model


V1.2.1
======

- Reverts to substitution method for partial fraction analysis

- Fixes factor_const and term_const


V1.2
====

- Add `discretize()` method for `TimeDomainExpression`

- Ignores `UnitStep` and conditional for Z-transform

- Scales `bilinear_transform()` by `1  / dt`

- Allows transformations from continuous-time to discrete-time

- Supports color arg for lollipop plots

- Fixes assumptions when scaling by a constant

- Adds Simpson, Euler, impulse-invariance, and matched-Z methods for discretization

- Generalizes `simplify_sin_cos`

- Adds include and includefile options for schtex

- Specifies voltage dir for Circuitikz

- Adds approximations for `exp`, `sinh`, `cosh`, `tanh`

- Fixes loop and nodal analysis in Laplace domain

- Improves simplification with complex conjugates

- Supports A and G two-ports for netlists

- Converts Greek names to symbols for schematics

- Adds `re` and `im` functions

- Speeds up inverse Laplace transform by computing residues by equating coefficients


V1.1
====

- Adds `loop_analysis` and `nodal_analysis` methods to `Circuit`

- Fixes creating two-port from netlist

- Improves Laplace transforms for convolutions

- Adds `Min` and `Max` functions

- Adds `solve()` method to `Expr` to solve expression

- Adds `solve()` methods to `ExprDict`, `ExprTuple`, and `ExprList` to solve system of equations

- Supports `AppliedUndef` for `Function`

- Uses `warn()` function throughout


V1.00
=====

- Overhauls `TwoPort` and associated classes

- Adds schematic support for two-ports

- Adds `solve()` to `ExprList` and `ExprTuple`

- Adds `Derivative`, `Integral`, and `Piecewise` functions

- Adds drawing hints to `Network` objects

- Fixes anonymous component names

- Adds MNA stamps for two-ports

- Adds `annotate_node_voltages()`, `annotate_voltages()`, and `annotate_currents()` methods

- Speeds up some Laplace Transforms

- Fixes odd bugs

- Fixes compatability with SymPy-1.9


V0.99
=====

- Separates state-space generation from state-space representation

- Adds discrete-time state-space representation `DTStateSpace`

- Adds creation of state-space models from transfer functions

- Adds state-space balancing

- Adds state-space model reduction

- Adds many DFTs

- Checks if have series L and independent current source for state-space generation

- Makes `Piecewise`, `Ne`, `Lt`, `Le`, `Gt`, `Ge` Lcapy functions

- Generalizes model discretization

- Adds matrix classes for discrete-time domains

- Adds Nichols plots

- Fixes printing of Piecewise

- Makes `conjugate` a method and adds `conj` as an attribute

- Fixes `evalf()`

- Adds `a` and `b` attributes for denominator and numerator coefficients


V0.98
=====

- Adds numerical filtering to `DLTIFilter`

- Normalizes a0 to 1 by default for `DLTIFilter`

- Add `subs()` method to `DLTIFilter`

- Fixes `subs()` method for `ExprDict`

- Adds inverse bilinear transform

- Adds `fval` and `cval` attributes to `ExprDict`, `ExprList`, and `ExprTuple`

- Ensures rationals converted to floats for `evalf()`

- Renames `form` with `layout` for network drawing

- Clarifies reasons why MNA fails

- Adds misc. bug fixes


V0.97
=====

- Adds many more DFTs

- Uses bilinear transform as default approach for `response()` method

- Preserves node order for loop finding

- Fixes domains of sequence elements

- Adds assumptions attribute to sequences

- Uses better naming for dummy variables



V0.96
=====

- Fixes `floatrat()` and `ratfloat()` expression methods

- Improves conversion of floats to rationals for `expr()`

- Ensures `evalf()` uses floats


V0.95
=====

- `expr()` handles `F` and `Omega` expressions

- Adds quantities and domains to sequences

- Adds domain argument to `seq`

- Fixes DFT caching

- Fixes plotting of discrete frequency expressions

- Supports sequences for `latex()` function


V0.94
=====

- Fixes plots

- Adds `dbmin` argument for frequency plots

- Fixes DTFTs

- Makes Heaviside and rect functions consistent with sign function

- Adds simplifications for Heaviside and rect functions

- Adds discrete-time rect and sign functions

- Warns if domain symbols are overridden

- Allows symbol redefinition

- Improves Nyquist plots



V0.93
=====

- Improves plotting dB-phase

- Plots Dirac deltas

- Speeds up plotting of frequency domain responses

- Adds Nyquist plots

- Fixes phasor transforms

- Evaluates Integrals, Sums, etc. before plotting

- Makes `is_complex` more robust

- Adds `pairs` argument to `ZPK()` to combine complex conjugates

- Adds `pairs` argument to `poles`, `zeros` and `roots` to combine complex conjugates

- Adds many more DTFTs

- Adds normalised frequency (F) and normalised angular frequency (Omega) domains

- Adds IDTFTs

- Ensures `dt` and `df` are positive

- Ensures `N` positive in DFT

- Adds generalized transformer infrastructure

- Fixes `dB`

- Warns about truncated sequences


V0.92
=====

- Fixes plotting frequency response

- Adds `norm` argument for frequency response plots

- Determines limit if NaN returned for `evaluate()`

- Adds `coth()` and `acoth()` functions

- Ensures `n` and `k` are integers

- Fixes `UnitStep` and `UnitImpulse`

- Adds `parameterize_ZPK()`

- Adds tutorial on expression manipulation

- Improves pole-zero plots


V0.91
=====

- Simplifies residues for better partial fractions

- Renames `DTFilter` to `DLTIFilter`

- Adds `DifferenceEquation` class

- Speeds up z-transforms

- Fixes stem plots for negative powers of n

- Ensures integer xticks for stem plots

- Adds `var` argument to `coeffs()` method for expressions

- Merges state space tests

- Changes behaviour of z-transform and DFT for sequences; they now return sequences

- Adds `expr` attribute for sequences

- Moves documentation to https:\\lcapy.readthedocs.org

- Improves pretty printing of sequences

- Adds `zeroextend()` method for sequences

- Adds `>>` and `<<` operators for sequences

- Uses attributes `extent` and `origin` for sequences

- Remove tests for deprecated ubuntu-16.04


V0.90
=====

- Adds call notation to access element of `Sequence`

- Adds `as_array()` method for `Sequence`

- Modifies `evaluate()` method for `Sequence`

- Adds `DTFilter`

- Adds `evalf` method for container classes

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

- Adds normalized and unnormalized versions of sinc

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

- V0.73 improves printing of Voltage and Current, adds phasor attributes to Voltage and Current, fixes magnitude and phase for Phasor, fixes printing of Greek symbols, tidies canonical representation, wraps R, X, B, G attributes for Immittance, doc reorganisation.

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

- V0.50 changes phasors to have a default angular frequency of omega_0 instead of omega to avoid confusion with angular frequency in Fourier transforms, adds preliminary phasor plots, improves noise signal classes, improves the infrastructure, and fixes many bugs.

- V0.49 adds mechanical components, better parameterization, faster partial fraction expansion, improved Z transforms, IIR difference equations, and differential equations.

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

- V0.36 improves handling of complex conjugate poles

- V0.34 switched to using setuptools and pushed to https::pypi.org

- V0.33 reworks expression printing infrastructure

- V0.32.3 introduces state-space analysis.  The API is experimental and may change.

- V0.32.0 changes the naming of symbolic values.  Previously R1 was converted to R_1 before being converted into a SymPy symbol.  This behaviour was not obvious for symbol substitution.  Now the symbol names are converted on printing.

- V0.31.0 reworks schematic drawing.  The syntax for chips has changed since there are no explicit nodes in the netlist.

- V0.30.0 tweaks the syntax to perform transformations based on the argument, e.g., V(s) or V(t)

- V0.28.0 works with Sympy 1.2

- V0.26.0 adds noise analysis

- V0.25.1 adds time-domain analysis for circuits without reactive components

- From version 0.25.0, Lcapy performs more comprehensive circuit analysis using combinations of DC, AC, and Laplace analysis.  This added functionality has resulted in a slight change of syntax.  cct.R1.V no longer prints the s-domain expression but the decomposition of a signal into each of the transform domains.

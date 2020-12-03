Lcapy is a Python package for linear circuit analysis.  It uses SymPy
for symbolic mathematics.

[![Travis-CI](https://api.travis-ci.org/mph-/lcapy.svg?branch=master)](https://travis-ci.org/mph-/lcapy)
[![Binder](http://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/mph-/lcapy/master)

Lcapy can analyse circuits described with netlists or by series/parallel combinations of components.

Comprehensive documentation can be found at http://lcapy.elec.canterbury.ac.nz


Circuit analysis
----------------

The circuit is described using netlists, similar to SPICE, with
arbitrary node names (except for the ground node which is labelled 0).
The netlists can be loaded from a file or created at run-time.  For
example:

    >>> from lcapy import Circuit, s
    >>> cct = Circuit("""
    ... Vs 2 0 {5 * u(t)}
    ... Ra 2 1
    ... Rb 1 0
    ... """)

The circuit can then be interrogated to determine branch currents,
branch voltages, and node voltages (with respect to the ground node 0).

    >>> cct[1].v
    >>> cct[2].v
    >>> cct.Ra.i
    >>> cct.Ra.V(s)


One-port networks
-----------------

One-port networks can be created by series and parallel combinations
of other one-port networks.  The primitive one-port networks are the
following ideal components:

- V independent voltage source
- I independent current source
- R resistor
- C capacitor
- L inductor

These components are converted to s-domain models and so capacitor and
inductor components can be specified with initial voltage and
currents, respectively, to model transient responses.

The components have the following attributes:

- Zoc open-circuit impedance
- Ysc short-circuit admittance
- Voc open-circuit voltage
- Isc short-circuit current

The component values can be specified numerically or symbolically
using strings, for example,

    >>> from lcapy import Vdc, R, L, C, s, t
    >>> R1 = R('R_1') 
    >>> L1 = L('L_1')
    >>> a = Vdc(10) + R1 + L1

Here a is the name of the network formed with a 10 V DC voltage source in
series with R1 and L1.

The s-domain open circuit voltage across the network can be printed with:

    >>> a.V(s)
    10/s

The time domain open circuit voltage is given by:

    >>> a.V(t)
    10

The s-domain short circuit current through the network can be printed with:

    >>> a.Isc(s)
    10/(L_1*s**2 + R_1*s)

The time domain short circuit current is given by:

    >>> a.Isc(t)
    10/R_1


Two-port networks
-----------------

One-port networks can be combined to form two-port networks.  Methods
are provided to determine transfer responses between the ports.

Here's an example of creating a voltage divider (L section)

    >>> from lcapy import *
    >>> a = LSection(R('R_1'), R('R_2'))


Limitations
-----------

1. Non-linear components cannot be modelled (apart from a linearisation around a bias point).

2. High order systems can go crazy.

3. Some two-ports generate singular matrices.


Schematics
----------

LaTeX schematics can be generated using circuitikz from the netlist.
Additional drawing hints, such as direction and size are required.

    >>> from lcapy import Circuit
    >>> cct = Circuit("""
    ... P1 1 0; down
    ... R1 1 3; right
    ... L1 3 2; right
    ... C1 3 0_1; down
    ... P2 2 0_2; down
    ... W 0 0_1; right
    ... W 0_1 0_2; right""")
    >>> cct.draw(filename='pic.tex')

In this example, P denotes a port (open-circuit) and W denotes a wire
(short-circuit).  The drawing hints are separated from the netlist
arguments by a semicolon.  They are a comma separated list of
key-value pairs except for directions where the dir keyword is
optional.  The symbol label can be changed using the l keyword; the
voltage and current labels are specified with the v and i keywords.
For example,

    >>> from lcapy import Circuit
    >>> cct = Circuit("""
    ... V1 1 0; down
    ... R1 1 2; left=2, i=I_1, v=V_{R_1}
    ... R2 1 3; right=2, i=I_2, v=V_{R_2}
    ... L1 2 0_1; down, i=I_1, v=V_{L_1}
    ... L2 3 0_3; down, i=I_1, v=V_{L_2}
    ... W 0 0_3; right
    ... W 0 0_1; left""")
    >>> cct.draw(scale=3, filename='pic2.svg')

The drawing direction is with respect to the positive node; i.e., the
drawing is performed from the positive to the negative node.  Since
lower voltages are usually lower in a schematic, then the direction of
voltage sources and ports is usually down.

By default, component (and current) labels are drawn above horizontal
components and to the right of vertical components.  Voltage labels
are drawn below horizontal components and to the left of vertical
components.

Node names containing a dot or underscore are not displayed.


Jupyter notebooks
-----------------

Lcapy can be used with [Jupyter Notebooks](https://jupyter.org/).  For a number of examples see https://github.com/mph-/lcapy/tree/master/doc/examples/notebooks .  These include:

- [AC analysis of a first-order RC filter](https://github.com/mph-/lcapy/blob/master/doc/examples/notebooks/RC-lpf1.ipynb)

- [A demonstration of the principle of superposition](https://github.com/mph-/lcapy/blob/master/doc/examples/notebooks/superposition2.ipynb)

- [Non-inverting operational amplifier](https://github.com/mph-/lcapy/blob/master/doc/examples/notebooks/opamp-noninverting-amplifier1.ipynb)

- [State-space analysis](https://github.com/mph-/lcapy/blob/master/doc/examples/notebooks/state-space1.ipynb)


Documentation
-------------

For comprehensive documentation, see http://lcapy.elec.canterbury.ac.nz

Alternatively, the documentation can be viewed in a web browser after
running 'make html' in the doc directory.

For another view on Lcapy see https://blog.ouseful.info/2018/08/07/an-easier-approach-to-electrical-circuit-diagram-generation-lcapy/


Testing
-------

The testsuite can be run using

    $ nosetests3 --pdb

Better still, use the --pdb option to enter the Python debugger on a failure:

    $ nosetests3 --pdb

To check for coverage use:

    $ nosetests3 -with-coverage --cover-package=lcapy --cover-html

and then view cover/index.html in a web browser.


Updates
-------

- Version 0.72 uses CI for docs plus many assorted bug fixes.

- Version 0.71 uses much faster matrix inversion (if sympy-1.8 installed) otherwise falls back on ADJ method  instead of the GE method which has a serious time regression with sympy-1.6.2

- Version 0.70 adds improved nodal and mesh analysis.

- Version 0.69 adds common-mode gain for opamps and polyphase-twoports.

- Version 0.67 adds time-stepping simulation, supernode detection, and polyphase circuits.

- Version 0.66 tidies up two-port parameters.  S and T parameters are
  added.  The A, B, G, H, Y, Z parameters are renamed to Aparams, etc. to avoid conflict with
  matrix transpose and Hermitian transpose operators.  issymmetrical, isshunt renamed to is_symmetrical,
  is+shunt, etc.   Eq, MatMul, MatAdd, Mul, Add functions added.  Expr.__getattr_ converts lists to ExprList.
  Adds symbols attribute to Matrix.  Ensures symbols in immitance default to complex.

- Version 0.65 introduces random networks.  Adds simplification for DiracDelta and Heaviside.  Adds node checking for Netlist methods.

- Version 0.64 adds wye-delta, delta-wye transformations.  Adds resistive companion models.  Fix state-space if have no sources.  Fixes assumption merging.  Adds verbatim argument for laplace_transform.   Simplifies mutual inductance.

- Version 0.63 fixes mirroring of opamps in schematics and introduces mirrorinputs option

- Version 0.62 adds search, save, annotate_voltage, annotate_current, kill_zero methods.  Fixes solve.

- Version 0.61 improves Laplace and z-transforms.

- Version 0.60 replaces DiracDelta with UnitImpulse and Heaviside with UnitStep for discrete-time expressions.

- Version 0.52 improves the component positioning algorithm for schematics.

- Version 0.51 improves the domain transformation infrastructure,

- Version 0.50 changes phasors to have a default angular frequeny of omega_0 instead of omega to avoid confusion with angular frequency in Fourier transforms, adds preliminary phasor plots, improves noise signal classes, improves the infrastructure, and fixes many bugs.

- Version 0.49 adds mechanical components, better parameterization, faster partical fraction expansion, improved Z transforms, IIR difference equations, and differential equations.

- Version 0.48 fixes z-transforms, adds better caching for Laplace and z-transforms, convert rational numbers to floats on schematics, fixes expr rpow.

- Version 0.47 introduces subs method for netlists, initialize method of netlists, better clarification for external programs, removes Y and Z methods for Circuits, removes anon ids from circuit components, adds remove_condition, force_causal, is_conditional, is_rational_function, is_strictly_proper, adds isoamp, inamp, and bug fixes

- Version 0.42 bug fixes for discrete-time signals

- Version 0.41 introduces experimental discrete-time signals

- Version 0.40 fixes schematics

- Version 0.39 miscellaneous bug fixes

- Version 0.38 reverts the experimental behaviour of 0.37.  Instead it introduces new classes for general immitances that tries to display them in the most suitable format.

- Version 0.37 changes the API for admittances and impedances.  The
  attributes Y and Z return the impedance in terms of omega rather
  than s as in the previous versions.  The old behaviour is provided
  with the Ys and Zs attributes (generalized admittance and
  impedance).  It also has better distinction between the impedance of
  a component and the driving point impedance.

- Version 0.36 has improved handling of complex conjugate poles.        

- Version 0.34 switched to using setuptools and pushed to https::pypi.org

- Version 0.33 reworked expression printing infrastructure

- Version 0.32.3 introduces state-space analysis.  The API is experimental and may change.

- Version 0.32.0 changes the naming of symbolic values.  Previously R1 was converted to R_1 before being converted into a SymPy symbol.  This behaviour was not obvious for symbol substitution.  Now the symbol names are converted on printing.

- Version 0.31.0 reworks schematic drawing.  The syntax for chips has changed since there are no explicit nodes in the netlist.

- Version 0.30.0 tweaks the syntax to perform transformations based on the argument, e.g., V(s) or V(t)

- Version 0.28.0 works with Sympy 1.2.

- Version 0.26.0 adds noise analysis.

- Version 0.25.1 adds time-domain analysis for circuits without reactive
components.

- From version 0.25.0, Lcapy performs more comprehensive circuit
analysis using combinations of DC, AC, and Laplace analysis.  This
added functionality has resulted in a slight change of syntax.
cct.R1.V no longer prints the s-domain expression but the
decomposition of a signal into each of the transform domains.


Copyright 2014--2020 Michael Hayes, UCECE

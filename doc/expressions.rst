.. _expressions:

===========
Expressions
===========

Lcapy expressions are similar to SymPy expressions except they have a
specific domain depending on the predefined domain variables `t`, `s`, `f`,
`omega`, and `jomega`.


Symbols
=======

Lcapy has a number of pre-defined constants, variables, and functions.


Constants
---------

- `pi` 3.141592653589793...

- `j`  :math:`\sqrt{-1}`

- `oo` infinity

- `zoo` complex infinity

Lcapy represents floating point numbers as rationals.  This ensures expected simplifications of expressions.    

.. _domainvariables:

Domain variables
----------------

Lcapy has five predefined domain variables:

- `s` Laplace domain complex frequency

- `f` Fourier domain frequency    

- `t` time
  
- `omega` Fourier domain angular frequency

- `jomega` Fourier domain angular frequency times `j`

A time-domain expression is produced using the `t` variable, for example,
  
   >>> v = exp(-3 * t) * u(t)

Similarly, a Laplace-domain expression is produced using the `s`
variable, for example,
  
   >>> V = s / (s**2 + 2 * s + 3)

The `discretetime` module adds additional domain variables `n`, `k`, and `z`, see :ref:`discrete_time`.
   
   
User defined symbols
--------------------

Symbols can also be created with Lcapy's `symbol` function:

   >>> tau = symbol('tau', real=True)

They are also implicitly created using Lcapy's `expr` function:
   
   >>> v = expr('exp(-t / tau) * u(t)')

Notes:

1. Symbols created with `symbol` and `expr` are assumed to be
   positive, unless explicitly specified not be.

2. Redefining a symbol does change the assumptions.  Instead, the symbol needs to be deleted with `symbol_delete` before being redefined.

3. There are restrictions on symbol names that can be used.  Currently, this excludes names that are Python keywords.  For example, `Is` is not allowed but `I_s` is valid.


.. _expressionsfunctions:
   
Mathematical functions
----------------------

Lcapy has the following built-in functions: `sin`, `cos`, `tan`,
`atan`, `atan2`, `gcd`, `exp`, `sqrt`, `log`, `log10`, `Heaviside`,
`H`, `u`, `DiracDelta`, `delta`, and `conjugate`.

Other SymPy functions can be converted to Lcapy functions using the
`Function` class, for example:

   >>> gamma = Function(sym.gamma)   


.. _expressionsrationalfunctions:
   
Rational functions
==================

Linear time-invariant systems have transfer functions that are rational functions; the ratio of two polynomials:

.. math::
   H(s) = \frac{N(s)}{D(s)},

The numerator can be found using the `N` attribute and denominator can
be found using the `D` attribute.
   

.. _expressionsresponses:
   
Responses
=========

Often, s-domain responses are rational functions but not always.  In general, Lcapy tries to interpret responses as

.. math::
   Y(s) = \sum_{i} \frac{N_i(s)}{D(s)} \exp(-s \tau_i),

where :math:`\tau_i` are time delays.   This representation is returned by the `as_sum()` method.  Note, these expressions cannot be represent in ZPK form.  The `D` attribute returns :math:`D(s)` and the `N` attribute returns

.. math::
   N(s) = \sum_{i} N_i(s) \exp(-s \tau_i).


.. _expressionsattributes:     

Attributes
==========

All Lcapy expressions have the following attributes (see :ref:`expressionsrationalfunctions` and :ref:`expressionsresponses` for definitions of numerator and denominator):

- `abs` returns absolute value

- `angle` returns phase angle (radians)
  
- `cartesian` returns expression in form `real + j * imag`

- `conjugate` returns complex conjugate

- `dB` returns magnitude in decibels: `20 * log10(magnitude)`

- `D` returns denominator

- `Ddegree` returns degree of denominator

- `denominator` returns denominator

- `degree` returns degree (order) of rational function (maximum of numerator and denominator degrees)
  
- `domain_label` returns string describing domain of expression
  
- `expr` returns the underlying SymPy expression
  
- `imag` returns imaginary part

- `is_ac` returns True if AC signal

- `is_causal` returns True if signal is causal, i.e, is 0 for :math:`t < 0`

- `is_conditional` returns True if expression is conditional, e.g., :math:`\exp(-t)\;\; \ge 0`

- `is_constant` returns True if expression constant

- `is_dc` returns True if DC signal    

- `is_number` returns True if expression is a number

- `is_rational_function` returns True if expression is a rational function

- `is_strictly_proper` returns True if degree of denominator greater than degree of numerator
  
- `label` returns string describing expression to use as a plot label

- `magnitude` returns absolute value  

- `N` returns numerator

- `Ndegree` returns degree of numerator    

- `numerator` returns numerator

- `phase` returns phase (radians)

- `phase_degrees` returns phase (degrees)    

- `polar` returns expression in form `mag * exp(j * phase)`

- `real` returns real part  

- `real_imag` returns expression in form `real + j * imag`

- `sign` returns sign

- `symbols` returns dictionary of symbols used in the expression keyed by their names

- `val` returns floating point number if expression can be evaluated

- `var` returns the underlying SymPy symbol representing the domain
    

.. _expressionsmethods:  

Methods
=======

Poles and zeros
---------------

- `coeffs()` returns list of coefficients if expression is a polynomial; the highest powers come first.  If the expression is a rational function use `.N.coeffs` or `.D.coeffs` for the numerator or denominator coefficients.

- `normcoeffs()` returns list of coefficients if expression is a polynomial; the highest powers come first.  The coefficients are normalised so the highest order coefficient is 1.  If the expression is a rational function use `.N.coeffs` or `.D.coeffs` for the numerator or denominator coefficients.

- `poles()` returns poles of expression as a dictionary or a list if the `aslist` argument is True.  Note, this does not always find all the poles.   

- `roots(s)` returns roots of expression as a dictionary or a list if the `aslist` argument is True.  Note, this does not always find all the roots.

- `zeros()` returns zeros of expression as a dictionary or a list if the `aslist` argument is True.  Note, this does not always find all the zeros.   
  

Miscellaneous
-------------

- `as_sum()` rewrite expression as a sum of terms where the denominator of each term has a common polynomial expression (see :ref:`expressionsresponses`).

- `divide_top_and_bottom(expr)` divides numerator and denominator by `expr`.

- `evalf()` returns floating point number if expression can be evaluated.
    
- `initial_value()` returns result at :math:`t = 0`.

- `factor_const()` factor into constant part and the rest.

- `factor_term()` split into constant part and the rest.    
  
- `final_value()` returns result at :math:`t = oo`.

- `multiply_top_and_bottom(expr)` multiplies numerator and denominator by `expr`.

- `rationalize_denominator` multiplies numerator and denominator by complex conjugate of denominator.

- `replace(query, value)` replace `query` with `value`.

  
.. _expressionsprinting:  
  
Formatting methods
------------------

Lcapy expressions can be displayed in many forms.  For example,
consider the s-domain rational-function:

   >>> H = 5 * (s**2 + 1) / (s**2 + 5*s + 4)     

The canonical form has a unity coefficient for the highest power in the denominator.  It is sometimes called polynomial form.
   
   >>> H.canonical()
     ⎛   2    ⎞ 
     ⎝5⋅s  + 5⎠   
   ────────────
    2          
   s  + 5⋅s + 4

There is a variation of the canonical form which has a unity coefficient for the highest power in the denominator and with constants factored in the numerator.   It is sometimes called gain-polynomial form.
   
   >>> H.canonical(factor_const=True)
      ⎛ 2    ⎞ 
    5⋅⎝s  + 1⎠ 
   ────────────
    2          
   s  + 5⋅s + 4

The general form of a rational function is shown as the ratio of two polynomials.   Unlike the canonical form, the coefficient for the highest power in the denominator may not be unity.
   
   >>> H.general()
        2      
     5⋅s  + 5  
   ────────────
    2          
   s  + 5⋅s + 4

The factored form show both the numerator and denominator polynomials  factored.  It is an alias for `ZPK` (zero-pole-gain) form.
   
   
   >>> H.factored()
   5⋅(s - ⅉ)⋅(s + ⅉ)
   ─────────────────
    (s + 1)⋅(s + 4) 

The partial fraction form has terms that are strictly proper.
    
   >>> H.partfrac()
           85          10   
   5 - ───────── + ─────────
       3⋅(s + 4)   3⋅(s + 1)

The `recippartfrac()` method generates a partial fraction expansion using the reciprocal of the variable:

   >>> H.recipartfrac()
   5       10          85    
   ─ - ───────── + ──────────
   4     ⎛    1⎞      ⎛1   1⎞
       3⋅⎜1 + ─⎟   48⋅⎜─ + ─⎟
         ⎝    s⎠      ⎝4   s⎠

       
The standard form expresses the rational function as the sum of a polynomial and a strictly proper rational function.
       
   >>> H.standard()
      25⋅s + 15      
   - ──────────── + 5
      2              
     s  + 5⋅s + 4    

The time constant form factors the rational function into gain-time-constant form.
   
   >>> H.timeconst()
   5⋅(-ⅉ⋅s + 1)⋅(ⅉ⋅s + 1)
   ──────────────────────
       ⎛s    ⎞           
     4⋅⎜─ + 1⎟⋅(s + 1)   
       ⎝4    ⎠           

The expanded canonical form expresses the rational function into the sum of rational functions where the numerator of each term contains a unique monomial.
       
   >>> H.expandcanonical()  
          2                   
       5⋅s             5      
   ──────────── + ────────────
    2              2          
   s  + 5⋅s + 4   s  + 5⋅s + 4


The `partfrac()` and `recippartfrac()` methods have a `combine_conjugates` argument.  If this is True, quadractic factors will not be split into two terms.  For example,

   >>> H = 5 / (s * (s**2 + 1))
   >>> H.partfrac()
         5           5       5
   - ───────── - ───────── + ─
     2⋅(s + ⅉ)   2⋅(s - ⅉ)   s
   >>> H.partfrac(combine_conjugates=True)
         5⋅s     5
      - ────── + ─
         2       s
        s  + 1    
  

Printing methods
----------------

- `pprint()` pretty print an expression

- `latex()`  convert an expression to LaTeX string representation

- `pretty()` convert an expression to a string with a prettified form

- `plot()` plot the expression, provided there are no free symbols
  

SymPy methods
-------------

If Lcapy does not have a method defined but the underlying SymPy
expression does, the SymPy method is used.  For example,

- `diff()`

- `simplify()`
  
   
Utility functions
=================

- `symbol()`  create a symbol

- `expr()` create an expression.  This can also create lists, tuples, and dictionaries of expressions.

Note, SymPy does not allow symbol names that are Python keywords.  For example,
`expr('is(t)')` fails.  A workaround is to use an underscore in the name, for example, `expr('i_s(t)')`.

- `simplify_terms()` expand expression into terms and simplify each term.

- `simplify_factor()` factor expression and simplify each factor.

- `limit()` compute a limit.  
  
Transformation and substitution
===============================      

Substitution and transformation use a similar syntax `V(arg)`.  If
`arg` is a domain variable `t`, `f`, `s`, `omega`, or `jomega`,
transformation is performed, otherwise substitution is performed.
This behaviour can be explicitly controlled using the `subs` and
`transform` methods, for example,

   >>> from lcapy import *
   >>> V1 = Voltage('3 * exp(-2 * t)')
   >>> V1.transform(s)
     3  
   ─────
   s + 2
   >>> V1.transform(t)
      -2⋅t
   3⋅e    
   >>> V1.subs(2)
      -4
   3⋅e  


Transformation
--------------


- `V(t)` returns the time domain transformation

- `V(f)` returns the Fourier domain transformation      

- `V(s)` returns the Laplace domain (s-domain) transformation

- `V(omega)` returns the angular Fourier domain transformation

- `V(jomega)` returns the angular Fourier domain transformation
  obtained from the Laplace domain transformation with :math:`s = j
  \omega`.

For example:

   >>> from lcapy import *
   >>> V1 = Voltage('3 * exp(-2 * t)')
   >>> V1(t)
      -2⋅t
   3⋅e    
   >>> V1(s)    
     3  
   ─────
   s + 2

  
Substitution
------------

Substitution replaces sub-expressions with new sub-expressions in an
expression.  It is most commonly used to replace the underlying
variable with a constant, for example,

   >>> a = 3 * s
   >>> b = a(2)
   >>> b
   6


Evaluation
----------
    
Evaluation is similar to substitution but requires all symbols in an
expression to be substituted with values.  The result is a numerical
answer.  The evaluation method is useful for plotting results.  For
example,

   >>> a = expr('t**2 + 2 * t + 1')
   >>> a.evaluate(0)
   1.0

The argument to `evaluate` can be a scalar, a tuple, a list, or a
NumPy array.  For example,

   >>> a = expr('t**2 + 2 * t + 1')
   >>> tv = np.linspace(0, 1, 5)
   >>> a.evaluate(tv)
   array([1.    , 1.5625, 2.25  , 3.0625, 4.    ])


Phasors
=======

Phasors represent signals of the form :math:`v(t) = A \cos(\omega t +
\phi)` as a complex amplitude :math:`X = A \exp(\mathrm{j} \phi)` where
:math:`A` is the amplitude, :math:`\phi` is the phase, and the angular
frequency, :math:`\omega`, is implied.

The signal :math:`v(t) = A \sin(\omega t)` has a phase
:math:`\phi=-\pi/2`.
      

.. _immitances:
      
Immitances
==========

Immitances (impedances and admittances) are represented using the
`Impedance` and `Admittance` classes.  They are primarily for internal
use.

Immitances can be initialised using either `omega` -domain or
`s` -domain expressions, for example:

   >>> Z1 = Impedance(5 * s)
   >>> Z2 = Impedance(5 * j * omega)

The impedance can be converted to a specific domain using a domain variable
as an argument.  For example,

   >>> Z1(s)
   >>> Z1(omega)

The time-domain representation of the immitance is the inverse Laplace
transform of the s-domain immittance, for example:

   >>> Impedance(1 / s)(t)
   Heaviside(t)
   >>> Impedance(1)(t)
   δ(t)
   >>> Impedance(s)(t)
    (1)    
   δ    (t)

The common way for creating an `Immitance` uses the `Y` or `Z` attribute of a
`Oneport` component, for example:

   >>> C(3).Z
   -ⅉ 
   ───
   3⋅ω

   >>> C(3).Z(s)
    1 
   ───
   3⋅s
   >>> C(3).Y(s)
   3⋅s

   

Immitance attributes
--------------------

- `B` susceptance

- `G` conductance    
  
- `R` resistance

- `X` reactance
  
- `Y` admittance

- `Z` impedance

Impedance is related to resistance and reactance by
  
:math:`Z = R + \mathrm{j} X`

Admittance is related to conductance and susceptance by      

:math:`Y = G + \mathrm{j} B`
        
Since admittance is the reciprocal of impedance,

:math:`Y = \frac{1}{Z} = \frac{R}{R^2 + X^2} - \mathrm{j} \frac{X}{R^2 + X^2}`

Thus

:math:`G = \frac{R}{R^2 + X^2}`

and

:math:`B = \frac{-X}{R^2 + X^2}`      
      
      
Note, at DC, when :math:`X = 0`, then :math:`G = 1 / R` and is
infinite for :math:`R= 0`.  However, if Z is purely imaginary, i.e,
:math:`R = 0` then :math:`G = 0`, not infinity as might be expected.
  

Immitance methods
-----------------
  
- `oneport()` returns a `Oneport` object corresponding to the immitance.  This may be a `R`, `C`, `L`, `G`, `Y`, or `Z` object.


Voltages and currents
=====================

Voltages and currents are represented using the `Voltage` and
`Current` classes.  These classes have similar behaviour; they
represent an arbitrary voltage or current signal as a superposition of
DC, AC, and transient signals.

For example, the following expression is a superposition of a DC
component, an AC component, and a transient component:

   >>> V1 = Voltage('1 + 2 * cos(2 * pi * 3 * t) + 3 * u(t)')

The signal can be converted to another domain using a domain variable
as an argument:

- `V1(t)` returns the time domain expression
- `V1(f)` returns the Fourier domain expression with linear frequency
- `V1(s)` returns the Laplace domain expression
- `V1(omega)` returns the Fourier domain expression with angular frequency
- `V1(jomega)` returns the Fourier domain expression with angular frequency    

Here are some examples,

   >>> V1(t)
   2⋅cos(6⋅π⋅t) + 3⋅u(t) + 1
   >>> V1(s)
     ⎛ 2       2⎞
   6⋅⎝s  + 24⋅π ⎠
   ──────────────
     ⎛ 2       2⎞
   s⋅⎝s  + 36⋅π ⎠
   >>> V1(jomega)
        ⎛   2       2⎞ 
   -6⋅ⅉ⋅⎝- ω  + 24⋅π ⎠ 
   ────────────────────
       ⎛   2       2⎞  
     ω⋅⎝- ω  + 36⋅π ⎠  



Voltage and current attributes
------------------------------

- `dc` returns the DC component
- `ac` returns a dictionary of the AC components, keyed by the frequency
- `transient` returns the time-domain transient component
- `is_dc` returns True if a pure DC signal
- `is_ac` returns True if a pure AC signal
- `is_transient` returns True if a pure transient signal
- `has_dc` returns True if has a DC signal
- `has_ac` returns True if has an AC signal
- `has_transient` returns True if has a transient signal


Voltage and current methods
---------------------------

- `oneport()` returns a `Oneport` object corresponding to the immitance.  This may be a `V` or `I` object.


Assumptions
===========

SymPy relies on assumptions to help simplify expressions.  In
addition, Lcapy requires assumptions to help determine inverse Laplace
transforms.

There are several attributes for determining assumptions:

- `is_dc` -- constant

- `is_ac` -- sinusoidal

- `is_causal` -- zero for :math:`t < 0`

- `is_real` -- real

- `is_complex` -- complex

- `is_positive` -- positive

- `is_integer` -- integer
    
For example:
  
   >>> t.is_complex  
   False
   >>> s.is_complex
   True
  

Assumptions for symbols
-----------------------

The more specific assumptions are, the easier it is for SymPy to solve
an expression.  For example,

   >>> C_1 = symbol('C_1', positive=True)

is more appropriate for a capacitor value than

   >>> C_1 = symbol('C_1', complex=True)


Notes:

   1. By default, the `symbol` and `expr` functions assume `positive=True` unless `real=True` or `positive=False` are specified.
   2. SymPy considers variables of the same name but different assumptions to be different.  This can cause much confusion since the variables look identical when printed.  To avoid this problem, Lcapy creates a symbol cache for each circuit.  The assumptions associated with the symbol are from when it is created.


The list of explicit assumptions for an expression can be found from
the `assumptions` attribute.  For example,

   >>> a = 2 * t + 3
   >>> a.assumptions
   {'real': True}

The `assumptions0` attribute shows all the assumptions assumed by SymPy.   

      
Assumptions for inverse Laplace transform
-----------------------------------------

Lcapy uses the :math:`\mathcal{L}_{-}` unilateral Laplace transform
(see :ref:`laplace_transforms`).  This ignores the function for
:math:`t <0` and thus the unilateral inverse Laplace transform thus
cannot determine the result for :math:`t <0` unless it has additional
information.  This is provided using assumptions:

-  `causal` says the signal is zero for :math:`t < 0`.

-  `ac` says the signal is sinusoidal.

-  `dc` says the signal is constant.

-  `damped_sin` says to write response of a second-order system as a damped sinusoid.
   
For example,

   >>> H = 1 / (s + 2)
   >>> H(t)
   ⎧ -2⋅t           
   ⎨e      for t ≥ 0
   ⎩                
   >>> H(t, causal=True)
    -2⋅t             
   e    ⋅Heaviside(t)

   >>> h = cos(6 * pi * t)
   >>> H = h(s)
   >>> H
       s     
   ──────────
    2       2
   s  + 36⋅π 
   >>> H(t)
   {cos(6⋅π⋅t)  for t ≥ 0
   >>> H(t, ac=True)
   cos(6⋅π⋅t)


Domain classes
==============

Lcapy has many expression classes, one for each combination of domain
(time, Fourier, Laplace, etc) and expression type (voltage, current,
impedance, admittance, transfer function).  For example, to represent
Laplace domain entities there are the following classes:

- `sExpr` generic Laplace-domain expression

- `Vs` Laplace-domain voltage

- `Is` Laplace-domain current

- `Hs` Laplace-domain transfer function

- `Ys` Laplace-domain admittance

- `Zs` Laplace-domain impedance


.. _parameterization:

Parameterization
================

Lcapy can parameterize a number of first order, second order, and third order s-domain expressions.  For example, 

   >>> H1 = 3 / (s + 2)
   >>> H1p, defs = H1.parameterize()
   >>> H1p
     K  
   ─────
   α + s
   >>> defs                                                                    
   {K: 3, alpha: 2}

Here `defs` is a dictionary of the parameter definitions.
   
The original expression can be obtained by substituting the parameter definitions into the parameterized expression:

   >>> H1p.subs(defs)                                                           
     3  
   ─────
   s + 2

Here's a second order example:

   >>> H2 = 3 / (s**2 + 2*s + 4)
   >>> H2p, defs = H2.parameterize()
   >>> H2p
              K         
   ───────────────────
     2               2
   ω₀  + 2⋅ω₀⋅s⋅ζ + s 
 
   >>> defs
   {K: 3, omega_0: 2, zeta: 1/2}

Second order systems can be parameterized in many ways.  Here's another example:

   >>> H2p, defs = H2.parameterize(zeta=False)
   >>> H2p
               K           
   ───────────────────────
     2    2              2
   ω₁  + s  + 2⋅s⋅σ₁ + σ₁ 

   >>> defs
   {K: 3, omega_1: √3, sigma_1: 1}


.. _network-synthesis:
   
Network synthesis
=================

Lcapy has experimental support for a number of network synthesis.
This produces a network model from an s-domain impedance or admittance
expression.  There are many methods, some specifically for simple
network such as R-L networks, and more general methods including
Foster and Cauer synthesis.

    >>> Z = Impedance(4*s**2 + 3 * s + one / 6) / (s**2 + 2 * s / 3)
    >>> n = Z.network('cauerI')
    >>> n
    ((C(1) + R(2)) | C(3)) + R(4)
    >>> n.Z(s).canonical()
    
    :math:`\frac{4 s^{2} + 3 s + \frac{1}{6}}{s^{2} + \frac{2 s}{3}}`

    >>> n.draw(form='ladder')
          
Note, in this example `one` is used to avoid generating a floating point number `1 / 6`.

    
  
SymPy
=====

The underlying SymPy expression can be obtained using the `expr`
attribute of an Lcapy expression.  For example,

   >>> a = 2 * t + 3
   >>> a.expr
   2⋅t + 3

The methods of the SymPy expression can be accessed from the Lcapy expression, for example,

   >>> a.as_ordered_terms()
   [2⋅t, 3]

Another example is accessing the SymPy symbol assumptions:

   >>> t.assumptions0
   {'commutative': True,
    'complex': True,
    'hermitian': True,
    'imaginary': False,
    'real': True}
   
Lcapy represents floating point numbers as rationals.  This ensures expected simplifications of expressions.



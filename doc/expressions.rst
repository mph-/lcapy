===========
Expressions
===========

Lcapy expressions are similar to SymPy expressions except they have a
specific domain depending on the predefined variables `t`, `s`, `f`,
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


Variables
---------

Lcapy has five predefined variables:

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
  
Symbols can also be created with Lcapy's `symbol` function:

   >>> tau = symbol('tau', real=True)

They are also implicitly created using Lcapy's `expr` function:
   
   >>> v = expr('exp(-t / tau) * u(t)')

Note, symbols created with `symbol` and `expr` are assumed to be
positive, unless explicitly specified not be.

There are restrictions on symbol names that can be used.  Currently, this excludes names that are Python keywords.  For example, `Is` is not allowed.


Mathematical functions
----------------------

Lcapy has the following built-in functions: `sin`, `cos`, `tan`,
`atan`, `atan2`, `gcd`, `exp`, `sqrt`, `log`, `log10`, `Heaviside`,
`H`, `u`, `DiracDelta`, `delta`, and `conjugate`.


Attributes
==========

All Lcapy expressions have the following attributes:

- `abs` return absolute value

- `angle` return phase angle (radians)
  
- `cartesian` return expression in form `real + j * imag`

- `conjugate` return complex conjugate

- `dB` return magnitude in decibels: `20 * log10(magnitude)`

- `degree` return degree of rational function (maximum of numerator and denominator degrees)
  
- `D` return denominator

- `Ddegree` return degree of denominator

- `denominator` return denominator

- `degree` return degree (order).  If expression is a rational function the degree is the maximum degree of the numerator and denominator.

- `domain_label` return string describing domain of expression
  
- `evalf` return floating point number if expression can be evaluated

- `imag` return imaginary part

- `is_ac` return True if AC signal

- `is_causal` return True if signal is causal, i.e, is 0 for :math:`t < 0`

- `is_constant` return True if expression constant

- `is_dc` return True if DC signal    

- `is_number` return True if expression is a number

- `label` return string describing expression to use as a plot label

- `magnitude` return absolute value  

- `N` return numerator

- `Ndegree` return degree of numerator    

- `numerator` return numerator

- `phase` return phase (radians)

- `phase_degrees` return phase (degrees)    

- `polar` return expression in form `mag * exp(j * phase)`

- `real` return real part  

- `real_imag` return expression in form `real + j * imag`

- `sign` return sign

- `strictly_proper` return True if degree of denominator greater than degree of numerator
  
- `symbols` return dictionary of symbols used in the expression keyed by their names
  

Methods
=======

Poles and zeros
---------------

- `coeffs` return list of coefficients if expression is a polynomial; the highest powers come first.  If the expression is a rational function use `.N.coeffs` or `.D.coeffs` for the numerator or denominator coefficients.

- `normcoeffs` return list of coefficients if expression is a polynomial; the highest powers come first.  The coefficients are normalised so the highest order coefficient is 1.  If the expression is a rational function use `.N.coeffs` or `.D.coeffs` for the numerator or denominator coefficients.

- `poles` return poles of expression as a dictionary or a list if the `aslist` argument is True.  Note, this does not always find all the poles.   

- `roots` return roots of expression as a dictionary or a list if the `aslist` argument is True.  Note, this does not always find all the roots.

- `zeros` return zeros of expression as a dictionary or a list if the `aslist` argument is True.  Note, this does not always find all the zeros.   
  

Miscellaneous
-------------

- `initial_value` return result at :math:`t = 0`

- `final_value` return result at :math:`t = oo`  

  

Formatting methods
------------------

Lcapy expressions can be displayed in many forms.  For example,
consider the s-domain rational-function:

   >>> H = 5 * (s**2 + 1) / (s**2 + 5*s + 4)     

   >>> H.canonical()
     ⎛   2    ⎞ 
     ⎝5⋅s  + 5⎠   
   ────────────
    2          
   s  + 5⋅s + 4

This has a unity coefficient for the highest power in the denominator.  It is sometimes called polynomial form.

   >>> H.canonical(factor_const=True)
      ⎛ 2    ⎞ 
    5⋅⎝s  + 1⎠ 
   ────────────
    2          
   s  + 5⋅s + 4

This has a unity coefficient for the highest power in the denominator and with constants factored in the numerator.   It is sometimes called gain-polynomial form.

   >>> H.general()
        2      
     5⋅s  + 5  
   ────────────
    2          
   s  + 5⋅s + 4

This is the general form of a rational function shown as the ratio of two polynomials.   Unlike the canonical form, the coefficient for the highest power in the denominator may not be unity.
   
   >>> H.factored()
   5⋅(s - ⅉ)⋅(s + ⅉ)
   ─────────────────
    (s + 1)⋅(s + 4) 

Here both the numerator and denominator polynomials are factored.  It is an alias for `ZPK` (zero-pole-gain) form.

   >>> H.partfrac()
           85          10   
   5 - ───────── + ─────────
       3⋅(s + 4)   3⋅(s + 1)

This splits the rational function into partial fraction form.
       
   >>> H.standard()
      25⋅s + 15      
   - ──────────── + 5
      2              
     s  + 5⋅s + 4    

This expresses the rational function into the sum of a polynomial and a strictly proper rational function.
     
   >>> H.timeconst()
   5⋅(-ⅉ⋅s + 1)⋅(ⅉ⋅s + 1)
   ──────────────────────
       ⎛s    ⎞           
     4⋅⎜─ + 1⎟⋅(s + 1)   
       ⎝4    ⎠           

This expresses the rational function in gain-time constant form.
       
   >>> H.expandcanonical()  
          2                   
       5⋅s             5      
   ──────────── + ────────────
    2              2          
   s  + 5⋅s + 4   s  + 5⋅s + 4


Printing methods
----------------

- `pprint` pretty print an expression

- `latex`  convert an expression to LaTeX string representation

- `pretty` convert an expression to a string with a prettified form

- `plot` plot the expression, provided there are no free symbols
  

SymPy methods
-------------

If Lcapy does not have a method defined but the underlying SymPy
expression does, the SymPy method is used.  For example,

- `diff`

- `simplify`
  
   
Utility functions
=================

- `symbol`  create a symbol

- `expr` create an expression.  This can also create lists, tuples, and dictionaries of expressions.

Note, sympy does not allow symbol names that are Python keywords.  For example,
`expr('is(t)')` fails.  A workaround is to use an underscore in the name, for example, `expr('i_s(t)')`.
  
  
Transformation and substitution
===============================      

Substitution and transformation use a similar syntax `V(arg)`.  If
`arg` is `t`, `f`, `s`, `omega`, or `jomega`, transformation is
performed, otherwise substitution is performed.  This behaviour can be
explicitly controlled using the `subs` and `transform` methods, for
example,

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


- :math:`V(t)` returns the time domain transformation

- :math:`V(f)` returns the Fourier domain transformation      

- :math:`V(s)` returns the Laplace domain (s-domain) transformation

- :math:`V(omega)` returns the angular Fourier domain transformation

- :math:`V(jomega)` returns the angular Fourier domain transformation
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
frequency, :math:`\omega` is implied.

The signal :math:`v(t) = A \sin(\omega t)` has a phase
:math:`\phi=-\pi/2`.
      

Immitances
==========

Immitances are represented using the `Impedance` and `Admittance` classes.


Immitance attributes
--------------------

- `B` susceptance

- `G` conductance    
  
- `R` resistance

- `X` reactance
  
- `Y` admittance

- `Z` impedance


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

The signal can be coverted to another domain using:

- `V1(t)` returns the time domain expression
- `V1(s)` returns the Laplace domain expression
- `V1(omega)` returns the Fourier domain expression with angular frequency
- `V1(f)` returns the Fourier domain expression with linear frequency


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

The unilateral Laplace transform ignores the function for :math:`t <
0`.  The unilateral inverse Laplace transform thus cannot determine
the result for :math:`t <0` unless it has additional information.
This is provided using assumptions:

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

Another example is accessing the assumptions that SymPy considers:

   >>> t.assumptions0
   {'commutative': True,
    'complex': True,
    'hermitian': True,
    'imaginary': False,
    'real': True}

Note, every real symbol is also considered complex although with no
imaginary part.  The proper way to test assumptions is to use the
attributes `is_complex`, `is_real`, etc.  For example,

   >>> t.is_real
   True
   >>> t.is_complex
   False

There can be difficulties with symbol assumptions when working with
SymPy.  By default sympy creates symbols with few assumptions, for example,

   >>> from sympy import Symbol
   >>> R1 = Symbol('R')
   >>> R1.assumptions0
   {'commutative': True}


On the other hand, by default, Lcapy assumes that symbols are
positive.  For example,

   >>> from lcapy import symbol
   >>> R2 = symbol('R')
   >>> R2.assumptions0
   {'commutative': True,
   'complex': True,
   'hermitian': True,
   'imaginary': False,
   'negative': False,
   'nonnegative': True,
   'nonpositive': False,
   'nonzero': True,
   'positive': True,
   'real': True,
   'zero': False}


Since `R1` and `R2` have different assumptions, SymPy considers them different symbols even though they are both defined as `R`.
   

Lcapy represents floating point numbers as rationals.  This ensures expected simplifications of expressions.


Be careful with zero substitutions; in general it is best to evaluate
a limit at zero.  For example,

    >>> x = symbol('x')
    >>> (x * (s + 1 / x)).subs(x, 0)
    0

    >>> x = symbol('x')
    >>> (x * (s + 1 / x)).limit(x, 0)
    1
    
Another approach is expand the expression to avoid the division:

    >>> x = symbol('x')
    >>> (x * (s + 1 / x)).expand().subs(x, 0)
    1

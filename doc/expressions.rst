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
--------

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


Mathematical functions
----------------------

Lcapy has the following built-in functions: `sin`, `cos`, `tan`,
`atan`, `atan2`, `gcd`, `exp`, `sqrt`, `log`, `log10`, `Heaviside`,
`H`, `u`, `DiracDelta`, `delta`, and `conjugate`.


Printing functions
------------------

- `pprint` pretty print an expression

- `latex`  convert an expression to LaTeX string representation

- `pretty` convert an expression to a string with a prettified form


Utility functions
-----------------

- `symbol`  create a symbol

- `expr` create an expression

  
Transformation and substitution
===============================      

Substitution and transformation use a similar syntax `V(arg)`.  If
`arg` is `t`, `f`, `s`, `omega`, or `jomega`, transformation is
performed, otherwise substitution is performed.  This behaviour can be
explicitly controlled using the `subs` and `transform` methods, for
example,

   >>> from lcapy import *
   >>> V1 = Vsuper('3 * exp(-2 * t)')
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
   >>> V1 = Vsuper('3 * exp(-2 * t)')
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

-  `causal` this says the signal is zero for :math:`t < 0`.

-  `ac` this says the signal is sinusoidal.

-  `dc` this says the signal is constant.


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


Classes
=======

Lcapy uses myriads of classes, one for each combination of domain
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

   >>>t.is_real
   True
   >>>t.is_complex
   False



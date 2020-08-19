===============
Troubleshooting
===============


Common problems
===============


Symbol aliases
--------------

SymPy treats symbols with different assumptions as different symbols
even if they have the same name.  To reduce this confusion, Lcapy
assumes that symbol names are not aliased.  It achieves this by
maintaining a list of defined symbols for each circuit.  However, it
is unaware of symbols created by SymPy.


Symbol assumptions
------------------

There can be difficulties with symbol assumptions when working with
SymPy.  By default SymPy creates symbols with few assumptions, for example,

   >>> from sympy import Symbol
   >>> R1 = Symbol('R')
   >>> R1.assumptions0
   {'commutative': True}

On the other hand, by default, Lcapy assumes that symbols are
positive (unless explicitly defined otherwise).  For example,

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


Since `R1` and `R2` have different assumptions, SymPy considers them different symbols even though they are both defined with the same name `R`.

Note, every real symbol is also considered complex although with no
imaginary part.  The proper way to test assumptions is to use the
attributes `is_complex`, `is_real`, etc.  For example,

   >>> t.is_real
   True
   >>> t.is_complex
   False


Zero substitution
-----------------

Be careful with zero substitutions.  For example, consider
    >>> x = symbol('x')
    >>> (x * (s + 1 / x)).subs(x, 0)
    0

In general it is safer (but slower) to evaluate a limit at zero.  

    >>> x = symbol('x')
    >>> (x * (s + 1 / x)).limit(x, 0)
    1
    
Another approach is expand the expression to avoid the division:

    >>> x = symbol('x')
    >>> (x * (s + 1 / x)).expand().subs(x, 0)
    1


Computation speed
-----------------

Lcapy can be slow for large problems due to the computational
complexity of the algorithms.  If speed is important, it is better to
substitute symbolic values with numerical values.

The results from slow computations are cached to improve the speed.

Some SymPy operations can take an unexpectedly long time, for example, `limit()`.


Floating point values
---------------------

Lcapy approximates floating point values as rational numbers.   This helps when simplifying expressions.  However, the conversion is an approximation.  For example, consider

   >>> s + 2 / 3

This becomes:
   
       3333333333333333
   s + ────────────────
       5000000000000000

In this case, Python evaluates 2 / 3 as a floating point number which is then converted to a rational number.  Unfortunately, this is not quite the same as 2 / 3.   The approximation can be avoided by bypassing the conversion of 2 / 3 to 0.666666666666, say by using:

   >>> expr('s + 2 / 3')
   s + 2/3

Another approach is to use:
 
   >>> s + one * 2 / 3
   s + 2/3



Debugging
=========


schtex
------

If `schtex` crashes, rerun it with the `-pdb` option.  This will enter the Python debugger when an unhandled exception is raised.


pdb method
----------

The Python debugger (pdb) can be entered using the `pdb()` method for many Lcapy classes.   For example, the inverse Laplace transform can be debugged for the expression `1 / (s + 2)` using:

   >>> (1 / (s + 2)).pdb().ILT()


debug method
------------

Expressions have a `debug()` method that prints the representation of the expresison, including symbol assumptions.  For example,

   >>> (1 / (s + 'a')).debug()                                                 
   sExpr(Pow(Add(s: {'nonpositive': False, 'nonzero': False, 'composite': False, 'real': False, 'negative': False, 'even': False, 'odd': False, 'prime': False, 'positive': False, 'nonnegative': False, 'integer': False, 'commutative': True, 'rational': False, 'zero': False, 'irrational': False},
              a: {'nonpositive': False, 'extended_nonpositive': False, 'hermitian': True, 'extended_positive': True, 'real': True, 'imaginary': False, 'negative': False, 'extended_real': True, 'infinite': False, 'extended_negative': False, 'extended_nonnegative': True, 'positive': True, 'nonnegative': True, 'extended_nonzero': True, 'finite': True, 'commutative': True, 'zero': False, 'complex': True, 'nonzero': True})
,
          -1)



   

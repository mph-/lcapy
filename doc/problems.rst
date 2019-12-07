========
Problems
========


Symbol aliases
==============

SymPy treats symbols with different assumptions as different symbols
even if they have the same name.  To reduce this confusion, Lcapy
assumes that symbol names are not aliased.  It achieves this by
maintaining a list of defined symbols for each circuit.  However, it
is unaware of symbols created by SymPy.


Symbol assumptions
==================

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
=================

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
=================

Lcapy can be slow for large problems due to the computational
complexity of the algorithms.  If speed is important, it is better to
substitute symbolic values with numerical values.

The results from slow computations are cached to improve the speed.

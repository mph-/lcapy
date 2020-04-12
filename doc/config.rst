=============
Configuration
=============


Printing
========

There are three ways to print an expression:

1. string (mostly for use in the debugger)
2. unicode (for the python/Ipython command line)
3. latex (for jupyter notebooks and LaTeX files)

For each way, the output representation of an expression can be changed using the dictionaries `str_expr_map`, `pretty_expr_map`, and `latex_expr_map`.  For example, to print the Heaviside "function" as an `H` instead of the default `u`, for the LaTeX representation use:

>>> from lcapy.config import latex_expr_map
>>> latex_expr_map[sym.Heaviside] = 'H'

Handling the the imaginary operator `sym.I` is trickier since this is ingrained as `i` in SymPy.  By default, Lcapy uses `j` but this can be forced to be printed as `i`, see `lcapy/config.py`.


Parsing
=======

Sympy has many predefined symbols.  Lcapy ignores the ones listed in `config.excludes`.  For example, `N` represents the numerical evaluation function.  If you wish to have this symbol, remove it from the `excludes` list, for example,

>>> from lcapy.config import excludes
>>> excludes.remove('N')

Here's another example, that removes the falling factorial function `ff`,

>>> from lcapy.config import excludes
>>> excludes.append('ff')


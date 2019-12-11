"""This module contains configuration information.

Copyright 2019 Michael Hayes, UCECE

"""

# SymPy symbols to exclude.  It might be easier to add the ones we want...
exclude  = ('I', 'C', 'O', 'S', 'N', 'E', 'E1', 'Q', 'beta', 'gamma', 'zeta',
            'Le', 'Lt', 'Ge', 'Gt', 'Ci')

# Aliases for SymPy symbols
aliases = {'delta': 'DiracDelta', 'step': 'Heaviside', 'u': 'Heaviside',
           'j': 'I'}

# String replacements when printing as LaTeX.  For example, SymPy uses
# theta for Heaviside's step.
latex_string_map = {r'\theta\left': r'u\left'}

import sympy as sym
str_expr_map = {sym.I: 'j'}

latex_expr_map = {sym.I: '\mathrm{j}'}

# Hack to pretty print i as j
junicode = '\u2149'
from sympy.printing.pretty.pretty_symbology import atoms_table
atoms_table['ImaginaryUnit'] = junicode

pretty_expr_map = {sym.I: junicode}

functions = ('heaviside', 'diracdelta', 'conjugate', 'sqrt', 'exp',
             'log', 'log10', 'sin', 'cos', 'tan', 'cot' 'asin',
             'acos', 'atan', 'atan2', 'acot', 'sinh', 'cosh', 'tanh', 'asinh',
             'acosh', 'atanh', 'gcd', 'abs')

# Words to format in Roman font for LaTeX expressions. 
subscripts = ('in', 'out', 'ref', 'rms', 'load', 'source', 'avg',
              'mean', 'peak', 'pk', 'pk-pk', 'pp', 'min', 'max', 'src', 'bat',
              'cc', 'ee', 'dd', 'ss', 'ih', 'il', 'oh', 'ol',
              'typ', 'pkg', 'comp', 'step')

words = ('alpha', 'beta', 'gamma', 'delta', 'eta', 'zeta', 'theta',
         'iota', 'kappa', 'mu', 'nu', 'omicron', 'pi', 'rho', 'sigma', 'tau',
         'upsilon', 'omega')


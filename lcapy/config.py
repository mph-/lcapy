"""This module contains configuration information.

Copyright 2019--2021 Michael Hayes, UCECE

"""

# SymPy symbols to exclude.  It might be easier to add the ones we want...
excludes = ['I', 'C', 'O', 'S', 'N', 'E', 'E1', 'Q', 'beta', 'gamma', 'zeta',
            'Le', 'Lt', 'Ge', 'Gt', 'Ci']

# Aliases for SymPy symbols
aliases = {'delta': 'DiracDelta', 'step': 'Heaviside', 'u': 'Heaviside',
           'H': 'Heaviside', 'j': 'I'}

import sympy as sym
str_expr_map = {sym.I: 'j'}

latex_expr_map = {sym.I: '\mathrm{j}', sym.Heaviside: 'u'}

# Hack to pretty print i as j
junicode = '\u2149'
from sympy.printing.pretty.pretty_symbology import atoms_table
atoms_table['ImaginaryUnit'] = junicode

pretty_expr_map = {sym.I: junicode, sym.Heaviside: 'u'}

# Words to format in Roman font for LaTeX expressions. 

functions = ('heaviside', 'diracdelta', 'conjugate', 'sqrt', 'exp',
             'log', 'log10', 'sin', 'cos', 'tan', 'cot' 'asin',
             'acos', 'atan', 'atan2', 'acot', 'sinh', 'cosh', 'tanh', 'asinh',
             'acosh', 'atanh', 'gcd', 'abs', 'unitimpulse', 'arg', 'sign',
             'rect', 'sinc')

subscripts = ('in', 'out', 'ref', 'rms', 'load', 'source', 'avg',
              'mean', 'peak', 'pk', 'pk-pk', 'pp', 'min', 'max', 'src', 'bat',
              'cc', 'ee', 'dd', 'ss', 'ih', 'il', 'oh', 'ol',
              'typ', 'pkg', 'comp', 'step')

words = ('alpha', 'beta', 'gamma', 'delta', 'eta', 'zeta', 'theta',
         'iota', 'kappa', 'mu', 'nu', 'omicron', 'pi', 'rho', 'sigma', 'tau',
         'upsilon', 'omega')

# Can be 'GE', 'LU', 'ADJ', 'LDL', 'CH', 'DM-GE', 'DM-LU', 'DM-charpoly'
# Note, the DM methods require the git version of sympy otherwise
# the fallback method is used.
try:
    from sympy.polys.domainmatrix import DomainMatrix
    matrix_inverse_method = 'DM-charpoly'
except:
    matrix_inverse_method = 'ADJ'
    
matrix_inverse_fallback_method = 'ADJ'


# Allow 1 + s, etc.
loose_units = True
# Print units with expression.
show_units = False
# Print units in abbreviataed form.
abbreviate_units = True
# Print units in canonical form.
canonical_units = False
# This cannot be enabled without loose units (e.g., s + 1, or s * 'R' * 'C' + 1 would fail).
check_units = True

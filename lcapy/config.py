# SymPy symbols to exclude.
exclude  = ('C', 'O', 'S', 'N', 'E', 'E1', 'Q', 'beta', 'gamma', 'zeta')

# Aliases for SymPy symbols
aliases = {'delta': 'DiracDelta', 'step': 'Heaviside', 'u': 'Heaviside',
           'j': 'I'}

# String replacements when printing as LaTeX.  For example, SymPy uses
# theta for Heaviside's step.
latex_string_map = {r'\theta\left': r'u\left'}

import sympy as sym
print_expr_map = {sym.I: 'j'}

# Hack to pretty print i as j
from sympy.printing.pretty.pretty_symbology import atoms_table
atoms_table['ImaginaryUnit'] = '\u2149'

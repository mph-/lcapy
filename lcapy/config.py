"""This module contains configuration information.

Copyright 2019--2021 Michael Hayes, UCECE

"""

# SymPy symbols to exclude.  It might be easier to add the ones we want...
# 'I' is removed to avoid conflict with oneport current source and
# since `j` is used instead.
from sympy.printing.pretty.pretty_symbology import atoms_table
import sympy as sym
excludes = ['I', 'C', 'O', 'S', 'N', 'E', 'E1', 'Q', 'beta', 'gamma', 'zeta',
            'Le', 'Lt', 'Ge', 'Gt', 'Ci']

# Aliases for SymPy symbols
aliases = {'delta': 'DiracDelta', 'step': 'Heaviside', 'u': 'Heaviside',
           'H': 'Heaviside', 'j': 'I'}

str_expr_map = {sym.I: 'j'}

latex_expr_map = {sym.I: '\mathrm{j}', sym.Heaviside: 'u'}

# Hack to pretty print i as j
junicode = '\u2149'
atoms_table['ImaginaryUnit'] = junicode

pretty_expr_map = {sym.I: junicode, sym.Heaviside: 'u'}

# Words to format in Roman font for LaTeX expressions.

functions = ('heaviside', 'diracdelta', 'conjugate', 'sqrt', 'exp',
             'log', 'log10', 'sin', 'cos', 'tan', 'cot' 'asin',
             'acos', 'atan', 'atan2', 'acot', 'sinh', 'cosh', 'tanh',
             'asinh', 'acosh', 'atanh', 'gcd', 'abs', 'unitimpulse',
             'arg', 'sign', 'rect', 'sinc', 'sincn', 'sincu', 'trap',
             'tri', 'ramp', 'rampstep', 'dtrect', 'dtsign', 'psinc')

subscripts = ('in', 'out', 'ref', 'rms', 'load', 'source', 'avg',
              'mean', 'peak', 'pk', 'pk-pk', 'pp', 'min', 'max', 'src', 'bat',
              'cc', 'ee', 'dd', 'ss', 'ih', 'il', 'oh', 'ol',
              'typ', 'pkg', 'comp', 'step')

greek_letter_names = ('alpha', 'beta', 'gamma', 'delta', 'epislon',
                      'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda',
                      'mu', 'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau',
                      'upsilon', 'phi', 'chi', 'psi', 'omega')

words = greek_letter_names

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
# Print units in abbreviated form.
abbreviate_units = True
# Print units in canonical form.
canonical_units = False
# This cannot be enabled without loose units (e.g., s + 1, or s * 'R' * 'C' + 1 would fail).
check_units = True
# Set to 'none' to print s - p; None to print SymPy default form -p + s.
# Note, SymPy will convert s - -1 to 1 + s unless Add is used with evaluate=False
printing_order = None

# Use X11 colours so that they will work with graphziv and dot2tex.
colours = {'startnode': 'green', 'endnode': 'green',
           'assignednode': 'Orchid1', 'unassignednode': 'SkyBlue1',
           'fixednode': 'Yellow', 'fixededge': 'Red', 'stretchedge': 'black'}

# Definition of H(0).  With H(0) = 0.5 then sgn(0) = 0 as expected by SymPy and NumPy
heaviside_zero = 0.5
# This is the common convention.
unitstep_zero = 1

implicit_default = 'sground'
autoground_default = 'sground'

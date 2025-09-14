"""This module contains configuration information.

Copyright 2019--2025 Michael Hayes, UCECE

"""

# See rcdefaults.py

# SymPy symbols to exclude.  It might be easier to add the ones we want...
# 'I' is removed to avoid conflict with oneport current source and
# since 'j' is used instead.
from sympy.printing.pretty.pretty_symbology import atoms_table
import sympy as sym
excludes = ['I', 'C', 'O', 'S', 'N', 'E', 'E1', 'Q', 'beta', 'gamma', 'zeta',
            'Le', 'Lt', 'Ge', 'Gt', 'Ci']

# Aliases for SymPy symbols
aliases = {'delta': 'DiracDelta', 'step': 'Heaviside', 'u': 'Heaviside',
           'H': 'Heaviside', 'j': 'I'}

# These mappings are set in rcdefaults
str_expr_map = {}
latex_expr_map = {}
pretty_expr_map = {}

# Words to format in Roman font for LaTeX expressions.

functions = ('heaviside', 'diracdelta', 'conjugate', 'sqrt', 'exp',
             'log', 'log10', 'sin', 'cos', 'tan', 'cot' 'asin',
             'acos', 'atan', 'atan2', 'acot', 'sinh', 'cosh', 'tanh',
             'asinh', 'acosh', 'atanh', 'gcd', 'abs', 'unitimpulse',
             'arg', 'sign', 'rect', 'sinc', 'sincn', 'sincu', 'trap',
             'tri', 'ramp', 'rampstep', 'dtrect', 'dtsign', 'psinc',
             'besselj', 'bessely', 'besseli', 'besselk', 'hankel1', 'hankel2',
             'erf', 'erfc', 'sec', 'csc', 'dirac')

subscripts = ('in', 'out', 'ref', 'rms', 'load', 'source', 'avg',
              'mean', 'peak', 'pk', 'pk-pk', 'pp', 'min', 'max', 'src', 'bat',
              'cc', 'ee', 'dd', 'ss', 'ih', 'il', 'oh', 'ol',
              'typ', 'pkg', 'comp', 'step')

greek_letter_names = ('alpha', 'beta', 'gamma', 'delta', 'epislon',
                      'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda',
                      'mu', 'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau',
                      'upsilon', 'phi', 'chi', 'psi', 'omega')

initialisms = ('OP', 'TP')

words = greek_letter_names + initialisms

# Note, SymPy will convert s - -1 to 1 + s unless Add is used with evaluate=False
printing_order = None

# Use X11 colours so that they will work with graphziv and dot2tex.
colours = {'startnode': 'green', 'endnode': 'green',
           'assignednode': 'Orchid1', 'unassignednode': 'SkyBlue1',
           'fixednode': 'Yellow', 'fixededge': 'Red', 'stretchedge': 'black'}

rcparams_user_filename = '~/.lcapy/lcapyrc'
rcparams_local_filename = 'lcapyrc'

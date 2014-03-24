"""
Mcircuit is a Python library for linear circuit analysis.  It uses SymPy
for symbolic mathematics.
"""

from __future__ import absolute_import, print_function

__version__ = "0.1-git"

import sys
if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.6 or above is required for Mcircuit.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed

del sys

from .mcircuit import *
from .netlist import *

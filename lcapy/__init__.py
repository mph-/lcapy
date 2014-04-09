"""
Lcapy is a Python library for linear circuit analysis.  It uses SymPy
for symbolic mathematics.

It can analyse circuits described with netlists using modified nodal
analysis.  See lcapy.netlist

Alternatively, it can analyse networks and circuits formed by
combining one, two, and three port networks.

Copyright 2014 Michael Hayes, UCECE
"""

from __future__ import absolute_import, print_function

__version__ = "0.2-git"

import sys
if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.6 or above is required for Lcapy.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed

del sys

from .mcircuit import *
from .netlist import *

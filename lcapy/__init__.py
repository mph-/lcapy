"""
Lcapy is a Python library for linear circuit analysis.  It uses SymPy
for symbolic mathematics.

Lcapy can analyse circuits described with netlists using modified nodal
analysis.  See lcapy.netlist

Alternatively, Lcapy can analyse networks and circuits formed by
combining one, two, and three port networks.

Copyright 2014 Michael Hayes, UCECE
"""

from __future__ import absolute_import, print_function

__version__ = "0.6-git"

import sys
if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.6 or above is required for Lcapy.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed

del sys


# List of symbols that get imported with 'from lcapy import *'
__all__ = []

# Add the modules that get searched to allow 'from lcapy import V'
# rather then having to specify module, 'from lcapy.oneport import V'
from . import core
__all__.extend(core.__all__)
from .core import *

from . import oneport
__all__.extend(oneport.__all__)
from .oneport import *

from . import twoport
__all__.extend(twoport.__all__)
from .twoport import *

from . import threeport
__all__.extend(threeport.__all__)
from .threeport import *

from . import netlist
__all__.extend(netlist.__all__)
from .netlist import *

from . import schematic
__all__.extend(schematic.__all__)
from .schematic import *

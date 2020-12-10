"""This module defines discrete-time symbols.

Copyright 2020 Michael Hayes, UCECE

"""

from .sym import symsymbol

nsym = symsymbol('n', real=True)
ksym = symsymbol('k', real=True)
zsym = symsymbol('z', real=False)

dt = symsymbol('Delta_t', real=True)
df = symsymbol('Delta_f', real=True)

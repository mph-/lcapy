"""This module defines discrete-time symbols.

Copyright 2020---2021 Michael Hayes, UCECE

"""

from .sym import symsymbol

nsym = symsymbol('n', integer=True)
ksym = symsymbol('k', integer=True)
zsym = symsymbol('z', real=False)
Omegasym = symsymbol('Omega', real=True)
Fsym = symsymbol('F', real=True)

dt = symsymbol('Delta_t', real=True, positive=True)
df = symsymbol('Delta_f', real=True, positive=True)

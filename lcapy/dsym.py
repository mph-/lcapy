"""This module defines discrete-time symbols.

Copyright 2020---2021 Michael Hayes, UCECE

"""

from .sym import symsymbol, domainsymbol

nsym = domainsymbol('n', integer=True)
ksym = domainsymbol('k', integer=True)
zsym = domainsymbol('z', real=False)

dt = symsymbol('Delta_t', real=True, positive=True)
df = symsymbol('Delta_f', real=True, positive=True)

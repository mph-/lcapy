from lcapy import Circuit

cct = Circuit()
#cct.add('V1 1 0')
cct.add('L1 1 0 L1')
cct.add('L2 2 0 L1*a**2')
cct.add('K1 L1 L2 k')
cct.add('R1 2 3')
cct.add('R2 3 4')
cct.add('E1 4 0 3 0 A')

H = cct.transfer(1, 0, 4, 0)
H.pprint()

import sympy as sym
expr = H.expr
sym.pprint(sym.limit(expr, 'A', sym.oo))

from lcapy import *
import sympy as sym

s, T = sym.symbols('s T')

H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6) * sym.exp(-s * T)

pprint(H)
pprint(ZPK(H))
pprint(partfrac(H))
pprint(inverse_laplace(H))


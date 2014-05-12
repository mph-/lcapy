from lcapy import *
from sympy import Idc

p = (10 + 20 * I)

H = 1 / ((s + p) * (s + p.conjugate()))

pprint(partfrac(H))

pprint(inverse_laplace(H))




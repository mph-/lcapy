from lcapy import *

s = sExpr.s

H = 5 * (s + 5) / ((s + 3) * (s + 3))

pprint(ZPK(H))

pprint(partfrac(H))

pprint(inverse_laplace(H))




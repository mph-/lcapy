from lcapy import *

s = sExpr.s

G = 1 / (s**2 + 5 * s + 6)

pprint(G)
pprint(partfrac(G))


H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6)

pprint(H)
pprint(ZPK(H))
pprint(partfrac(H))


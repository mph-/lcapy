from mcircuit import *

s = sExpr.s

g = 1 / (s**2 + 5 * s + 6)

pprint(g)
pprint(partfrac(g))


h = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6)

pprint(h)
pprint(ZPK(h))
pprint(partfrac(h))


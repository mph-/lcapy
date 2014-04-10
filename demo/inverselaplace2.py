from lcapy import *

H = 5 * (s + 5) * (s - 4) / (s**2 + 5 * s + 6)

pprint(partfrac(H))

pprint(inverse_laplace(H))




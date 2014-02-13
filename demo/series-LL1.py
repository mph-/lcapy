from mcircuit import V, R, L, C, pprint

L1 = L(10)
L2 = L(20)

L3 = L1 + L2

pprint(L3)

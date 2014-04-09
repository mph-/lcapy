from lcapy import *

a = Chain(Shunt(R(10)), Series(L(5)))

a

print(a)
pprint(a)

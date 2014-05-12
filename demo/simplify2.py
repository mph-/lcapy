from lcapy import *

a = Shunt(R(10))

b = Shunt(R(10))

c = a.chain(b)

print(c)
print(c.simplify())


d = Shunt(R(10) + Vdc(10))

e = Shunt(R(5) + Vdc(5))

f = d.chain(e)

print(f)
print(f.simplify())



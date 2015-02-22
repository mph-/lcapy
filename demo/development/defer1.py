from lcapy import *

a = C(1)

b = L(3)

c = Par(a, b)

print(c)

print(c.simplify())



d = C(1)

e = C(3)

f = Par(d, e)

print(f)

print(f.simplify())


a = Vdc(1)

b = L(3)

c = Par(a, b)

print(c)

print(c.simplify())



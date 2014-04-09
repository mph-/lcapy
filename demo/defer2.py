from lcapy import *

a = C(1)

b = L(3)

c = Ser(a, b)

print(c)

print(c.simplify())



d = C(1)

e = C(3)

f = Ser(d, e)

print(f)

print(f.simplify())


a = V(1)

b = L(3)

c = Ser(a, b)

print(c)

print(c.simplify())


a = V(1)

b = V(3)

c = Ser(a, b)

print(c)

print(c.simplify())

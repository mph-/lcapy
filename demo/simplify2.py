from mcircuit import *

a = Shunt(R(10))

b = Shunt(R(10))

c = a.chain(b)

print(c)
print(c.simplify())


d = Shunt(R(10) + V(10))

e = Shunt(R(5) + V(5))

f = d.chain(e)

print(f)
print(f.simplify())



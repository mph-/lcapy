from mcircuit import *

a = C(1)

b = L(3)

c = Par(a, b)

print(c)

print(c.simplify())



a = C(1)

b = C(3)

c = Par(a, b)

print(c)

print(c.simplify())


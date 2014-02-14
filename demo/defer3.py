from mcircuit import *

a = C(1)
b = C(2)
c = C(3)
d = C(4)

f = Ser(c, d)

e = Par(Ser(a, b), Ser(c, d))

print(e)

print(e.simplify())

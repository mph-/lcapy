from lcapy import *

a = C(1)
b = C(2)
c = C(3)
d = C(4)

f = Ser(c, d)

e = Par(Ser(a, b), Ser(c, d))

print(e)

print(e.simplify())


a = Vdc(1)
b = C(2)
c = C(3)
d = C(4)

e = Par(Ser(a, b), Ser(c, d))

print(e)

print(e.simplify())



a = Vdc(1)
b = C(2)
c = C(3)
d = C(4)

e = Ser(Par(a, b), Par(c, d))

print(e)

print(e.simplify())



a = Vdc(1)
b = C(2)
c = C(3)
d = C(4)

e = Par(Par(a, b), Par(c, d))

print(e)

print(e.simplify())

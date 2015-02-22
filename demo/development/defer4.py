from lcapy import *

a = Vdc(1)
b = C(2)
c = Vdc(3)
d = C(4)

e = Par(Ser(a, b), Ser(c, d))

ee = e.simplify()

print(e)

print(ee)



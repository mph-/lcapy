from lcapy import *

a = V(1)
b = C(2)
c = V(3)
d = C(4)

e = Par(Ser(a, b), Ser(c, d))

ee = e.simplify()

print(e)

print(ee)



from mcircuit import *

a = TSection(R(10), R(20), R(30))
print(a)

b = a.Pisection()
print(b)

c = b.Tsection()
print(c)

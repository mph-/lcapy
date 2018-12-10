from lcapy import Vdc, R, L, C, I, pprint, Idc

L1 = L(10)
L2 = L(20)

print('Series inductors')
L3 = L1 + L2
pprint(L3)
pprint(L3.simplify())
pprint(L3.norton().simplify())


print('Parallel inductors')
L4 = L1 | L2
pprint(L4)
pprint(L4.simplify())
pprint(L4.norton().simplify())


C1 = C(10)
C2 = C(20)

C3 = C1 + C2

print('Series capacitors')
pprint(C3)
pprint(C3.simplify())
pprint(C3.norton().simplify())


C4 = C1 | C2

print('Parallel capacitors')
pprint(C4)
pprint(C4.simplify())
pprint(C4.norton().simplify())



V1 = Vdc(10)
V2 = Vdc(20)

V3 = V1 + V2

print('Series voltage sources')
pprint(V3)
pprint(V3.simplify())
pprint(V3.norton().simplify())
pprint(V3.thevenin().simplify())


I1 = Idc(10)
I2 = Idc(20)

I3 = I1 | I2

print('Parallel current sources')
pprint(I3)
pprint(I3.simplify())
pprint(I3.norton().simplify())
pprint(I3.thevenin().simplify())

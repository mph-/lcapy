from mcircuit import R, L, C, LSection, pprint, ZPK

C1 = C('C_1')
L1 = L('L_1')
C2 = C('C_2')
L2 = L('L_2')

a = LSection(C1 | L1, C2 | L2)

Av = a.Vtransfer

pprint(Av.canonical())

pprint(Av.ZPK())

pprint(Av.poles())

pprint(Av.zeros())


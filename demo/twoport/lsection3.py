from lcapy import R, L, C, LSection, pprint, ZPK

C1 = C('C_1')
L1 = L('L_1')
R1 = R('R_1')

a = LSection(R1, C1 | L1)

Av = a.Vtransfer

pprint(Av.canonical())

pprint(Av.ZPK())

pprint(Av.poles())

pprint(Av.zeros())


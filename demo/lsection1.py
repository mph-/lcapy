from mcircuit import R, L, C, LSection, pprint, ZPK

R1 = R('R_1')
R2 = R('R_2')

a = LSection(R1, R2)

Av = a.Vtransfer

pprint(Av.canonical)

pprint(Av.ZPK)

pprint(Av.poles())

pprint(Av.zeros())


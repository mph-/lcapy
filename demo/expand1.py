from mcircuit import R, L, C, Xtal, pprint

xtal = Xtal('C_0', 'R_1', 'L_1', 'C_1')
pprint(xtal)

foo = xtal | C('C_2')
pprint(foo)
pprint(foo.simplify())
pprint(foo.expand())
pprint(foo.expand().simplify())




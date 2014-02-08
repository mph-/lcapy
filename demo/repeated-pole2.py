from mcircuit import invLT
import sympy as sym

t, s = sym.symbols('t s')

H = (3 * s + 2) / ((s + 1)**2)

sym.pprint(sym.Eq(sym.sympify('H'), H))

sym.pprint(sym.Eq(sym.sympify('h'), invLT(H, s, t)))


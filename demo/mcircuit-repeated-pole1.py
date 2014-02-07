from mcircuit import invLT
import sympy as sym

t, s = sym.symbols('t s')

H = 2 / (s**3 + 12 * s**2 + 36 * s)

sym.pprint(sym.Eq(sym.sympify('H'), H))

sym.pprint(sym.Eq(sym.sympify('h'), invLT(H, s, t)))


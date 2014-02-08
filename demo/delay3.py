from mcircuit import invLT
import sympy as sym

t, s = sym.symbols('t s')

H = (sym.exp(10 * s) * (s + 3))/ (s**2 + 5 * s + 6)

sym.pprint(sym.Eq(sym.sympify('H'), H))

sym.pprint(sym.Eq(sym.sympify('h'), invLT(H, s, t)))


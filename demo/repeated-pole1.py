from lcapy import inverse_laplace
import sympy as sym

t, s = sym.symbols('t s')

H = 2 / (s**3 + 12 * s**2 + 36 * s)

sym.pprint(sym.Eq(sym.sympify('H'), H))

sym.pprint(sym.Eq(sym.sympify('h'), inverse_laplace(H, t, s)))


from lcapy import inverse_laplace, partfrac
import sympy as sym

t, s = sym.symbols('t s')

H = (s**2 + 4) / (s**2 + 3 * s + 2)

sym.pprint(sym.Eq(sym.sympify('H'), H))

sym.pprint(sym.Eq(sym.sympify('H'), partfrac(H, s)))

sym.pprint(sym.Eq(sym.sympify('h'), inverse_laplace(H, t, s)))


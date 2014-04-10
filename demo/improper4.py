from lcapy import inverse_laplace, partfrac
import sympy as sym

A, B, t, s = sym.symbols('A B t s')

H = (s**2 + A) / (s**2 + B * s + 2)

sym.pprint(sym.Eq(sym.sympify('H'), H))

sym.pprint(sym.Eq(sym.sympify('H'), partfrac(H, s)))

sym.pprint(sym.Eq(sym.sympify('h'), inverse_laplace(H, t, s)))


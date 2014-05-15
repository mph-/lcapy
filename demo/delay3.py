from lcapy import inverse_laplace,  s
import sympy as sym

H = (sym.exp(-10 * s) * (s + 4))/ (s**2 + 5 * s + 6)

sym.pprint(sym.Eq(sym.sympify('H'), H.expr))

sym.pprint(sym.Eq(sym.sympify('h'), H.inverse_laplace().expr))


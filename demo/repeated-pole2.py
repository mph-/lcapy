from mcircuit import inverse_laplace
import sympy as sym

t, s = sym.symbols('t s')

H = (3 * s + 2) / ((s + 1)**2)

sym.pprint(sym.Eq(sym.sympify('H'), H))

sym.pprint(sym.Eq(sym.sympify('h'), inverse_laplace(H, s, t)))


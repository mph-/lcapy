import sympy as sym


def residues(expr, var):
    numer, denom = expr.as_numer_denom()
    poles = sym.roots(sym.Poly(denom, var))
    return [sym.residue(expr, var, p) for p in poles]


def poles(expr, var):
    numer, denom = expr.as_numer_denom()
    poles = sym.roots(sym.Poly(denom, var))
    return poles

s = sym.symbols('s')

z = (s**3 + 5)/((s**4 - 1)*(s + 1))
print(poles(z, s))
print(residues(z, s))


from scipy.signal import residue
from scipy import poly1d

N = poly1d((50, 0, 5), variable='s')
D = poly1d((10, 2, 1, 0), variable='s')
R, P, Q = residue(N.c, D.c)
print P
print R


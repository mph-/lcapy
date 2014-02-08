import sympy as sym

def poles(expr, var):
    numer, denom = expr.as_numer_denom()
    poles = sym.roots(sym.Poly(denom, var))
    return poles


def residues(expr, var):

    P = poles(expr, var)
    F = []
    R = []
    for p in P:
        
        # Number of occurences of the pole.
        N = P[p]

        D = var - p

        if N == 1:
            tmp = expr * D
            F.append(D)
            R.append(sym.limit(tmp, var, p))
            continue

        # Handle repeated poles.
        expr2 = expr * D ** N
        for n in range(1, N + 1):
            M = N - n
            F.append(D ** n)
            R.append(sym.limit(sym.diff(expr2, var, M), var, p) / sym.factorial(M))

    return F, R


def partfrac(expr, var):

    F, R = residues(expr, var)

    expr = 0
    for f, r in zip(F, R):
        expr = expr + r / f

    return expr


s = sym.symbols('s')

z = (50*s**2 + 5)/(10*s**3 + 2*s**2 + s)
#print(residues(z, s))
sym.pprint(partfrac(z, s))

z = (3 * s + 2) / (s**2 + 2 * s + 1)
#print(residues(z, s))
sym.pprint(partfrac(z, s))

z = 2 / (s**3 + 12 * s**2 + 36 * s)
#print(residues(z, s))
sym.pprint(partfrac(z, s))




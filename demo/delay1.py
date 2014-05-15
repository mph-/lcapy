from lcapy import s, sExpr
import sympy as sym

H = sExpr(sym.exp(-10 * s))

sym.pprint(sym.Eq(sym.sympify('H'), H.expr))

sym.pprint(sym.Eq(sym.sympify('h'), H.inverse_laplace().expr))


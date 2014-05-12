from lcapy import *
import sympy as sym

a, b, c, d, t = sym.symbols('a b c d t', real=True)

r = a + sym.I * b
p = c + sym.I * d

h = r * sym.exp(p * t) + r.conjugate() * sym.exp(p.conjugate() * t)






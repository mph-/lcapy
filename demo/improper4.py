from lcapy import s
import sympy as sym

A, B = sym.symbols('A B')

H = (s**2 + A) / (s**2 + B * s + 2)

H.pprintans('H(s)')

H.partfrac().pprintans('H(s)')

H.inverse_laplace().pprintans('h(t)')








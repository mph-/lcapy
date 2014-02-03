from msignal.mcircuit import V, R, L, C
import sympy as sym

a = V(10) + R('R') + C('C') + L('L')

a.Isc.pprint()

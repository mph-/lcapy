from lcapy import V, R, L, C
import sympy as sym

R1 = R('R') 
L1 = L('L')

a = V(10) + R1 + L1

a.Isc.pprint()

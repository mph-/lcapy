from lcapy import Vstep, R, L, C
import sympy as sym

R1 = R('R') 
L1 = L('L')

a = Vstep(10) + R1 + L1

a.Isc.pprint()

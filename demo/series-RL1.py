from lcapy import Vdc, R, L, C
import sympy as sym

R1 = R('R') 
L1 = L('L')

a = Vdc(10) + R1 + L1

a.Isc.pprint()

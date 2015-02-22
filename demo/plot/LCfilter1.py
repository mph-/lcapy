from lcapy import Circuit, j, f, pi
from numpy import linspace

cct = Circuit()

cct.add('P1 1 0; down')
cct.add('L1 1 2 0.1; right') 
cct.add('C1 2 3 1e-3; right') 
cct.add('R1 3 4 10; down') 
cct.add('W1 0 4; right') 
H = cct.transfer(1, 0, 3, 4)
A = H(j * 2 * pi * f)

fv = linspace(0, 100, 400)

A.magnitude.dB.plot(fv)
A.phase.plot(fv)





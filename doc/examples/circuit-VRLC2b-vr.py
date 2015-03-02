from lcapy import Circuit
from matplotlib.pyplot import savefig, show
from numpy import linspace


cct = Circuit()

cct.add('V1 1 0 dc')
cct.add('L1 1 2')
cct.add('C1 2 3')
cct.add('R1 3 0')

Vr = cct.V[3]
vr = Vr.transient_response()

vr.pprintans('v_R(t)')

if False:
    tv = linspace(0, 2, 400)
    vr.plot(tv)
    
    show()
    savefig('circuit-VRLC2b-vr.png')

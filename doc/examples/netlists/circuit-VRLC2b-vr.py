from lcapy import Circuit, t
from matplotlib.pyplot import savefig
from numpy import linspace


cct = Circuit()

cct.add('V1 1 0 step')
cct.add('L1 1 2 L1 0')
cct.add('C1 2 3')
cct.add('R1 3 0')

Vr = cct[3].V
vr = Vr(t)

vr.pprintans('v_R(t)')

vr = vr.subs({'L1':1, 'C1':0.1, 'R1':10, 'V1':1})

if True:
    tv = linspace(0, 2, 400)
    vr.plot(tv)
    
    savefig('circuit-VRLC2b-vr.png')

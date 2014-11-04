from lcapy import Circuit

cct = Circuit()

cct.add('V1 1 0 dc')
cct.add('L1 1 2')
cct.add('C1 2 3')
cct.add('R1 3 0')

Vr = cct.V[3]
vr = Vr.transient_response()

vr.pprintans('v_R(t)')

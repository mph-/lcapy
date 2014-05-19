from lcapy import Circuit

cct = Circuit()

cct.net_add('V1 1 0 dc')
cct.net_add('L1 1 2')
cct.net_add('C1 2 3')
cct.net_add('R1 3 0')

Vr = cct.V[3]
vr = Vr.transient_response()

vr.pprintans('v_R(t)')

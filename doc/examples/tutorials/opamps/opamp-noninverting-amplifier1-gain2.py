from lcapy import Circuit, t, oo

a = Circuit('opamp-noninverting-amplifier1.sch')

Vo = a[1].V(t)
Vo.pprint()
Vo.limit('Ad', oo).pprint()

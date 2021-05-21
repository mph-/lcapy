from lcapy import Circuit, t, oo

a = Circuit('opamp-transimpedance-amplifier1.sch')

Vo = a[1].V(t)
Vo.pprint()


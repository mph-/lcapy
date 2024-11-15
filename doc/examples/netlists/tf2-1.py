from lcapy import Circuit

a = Circuit('tf2.sch')

H1 = a.transfer('V1', 'R3')

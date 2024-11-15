from lcapy import Circuit

a = Circuit('tf2.sch')

H2 = a.transfer('V2', 'R3')

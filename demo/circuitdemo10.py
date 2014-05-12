from lcapy import pprint, Circuit

cct = Circuit('V R C')

cct.net_add('Vs 1 0') 
cct.net_add('R1 1 2') 
cct.net_add('C1 2 0 C1 Vc') 

pprint(cct.V)

pprint(cct.I)

pprint(cct.Voc(2, 0))

pprint(cct.Isc(2, 0))




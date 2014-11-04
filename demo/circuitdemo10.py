from lcapy import pprint, Circuit

cct = Circuit('V R C')

cct.add('Vs 1 0 dc') 
cct.add('R1 1 2') 
cct.add('C1 2 0 C1 Vc') 

pprint(cct.V)

pprint(cct.I)

pprint(cct.Voc(2, 0))

pprint(cct.Isc(2, 0))




from lcapy import pprint, Circuit

cct = Circuit('V R C')

cct.add('Vs 1 0 dc') 
cct.add('R1 1 2') 
cct.add('C1 2 0 1 -5') 


pprint(cct.V)

pprint(cct.I)




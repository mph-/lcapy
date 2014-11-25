from lcapy import pprint, Circuit

cct = Circuit()

# Inverting opamp
cct.add('Vs 1 0 step') 
cct.add('R1 1 2') 
cct.add('R2 3 2') 
cct.add('E1 3 0 2 0 1e6') 



pprint(cct.V)

pprint(cct.I)




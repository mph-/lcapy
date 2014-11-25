from lcapy import pprint, Circuit

cct = Circuit()

cct.add('Vs 1 0 step 10') 
cct.add('R1 2 1') 
cct.add('R2 3 2') 
cct.add('E1 3 0 2 0 1e6') 



pprint(cct.V)

pprint(cct.I)




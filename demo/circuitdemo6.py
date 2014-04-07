from mcircuit import pprint, Circuit

cct = Circuit('Non-inverting opamp')

cct.net_add('Vs 1 0 10') 
cct.net_add('Ri 1 0') 
cct.net_add('R1 2 0') 
cct.net_add('R2 3 2') 
cct.net_add('E1 3 0 2 1 1e6') 

cct.analyse()

pprint(cct.V)

pprint(cct.I)




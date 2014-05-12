from lcapy import pprint, Circuit

cct = Circuit('V R C')

cct.net_add('Vs 1 0') 
cct.net_add('R1 1 2') 
cct.net_add('C1 2 0 1 -5') 


pprint(cct.V)

pprint(cct.I)




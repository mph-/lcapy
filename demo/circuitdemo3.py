from lcapy import pprint, Circuit

cct = Circuit('V R C')

cct.net_add('Vs 1 0 10') 
cct.net_add('R1 1 2 5') 
cct.net_add('C1 2 0 1') 
cct.analyse()

pprint(cct.V)

pprint(cct.I)




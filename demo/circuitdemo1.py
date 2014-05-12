from lcapy import pprint, Circuit

cct = Circuit('Voltage divider')

cct.net_add('Vs 1 0') 
cct.net_add('Ra 1 2') 
cct.net_add('Rb 2 0') 


pprint(cct.V)

pprint(cct.I)

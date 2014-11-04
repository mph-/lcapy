from lcapy import pprint, Circuit

cct = Circuit('Voltage divider')

cct.add('Vs 1 0 dc') 
cct.add('Ra 1 2') 
cct.add('Rb 2 0') 


pprint(cct.V)

pprint(cct.I)

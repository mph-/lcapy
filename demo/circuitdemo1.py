from lcapy import pprint, Circuit

cct = Circuit()

cct.add('Vs 1 0 step') 
cct.add('Ra 1 2') 
cct.add('Rb 2 0') 


pprint(cct.V)

pprint(cct.I)

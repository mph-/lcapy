from lcapy import pprint, Circuit

cct = Circuit('V R C')

cct.add('V1 1 0 ac 10 100') 
cct.add('R1 1 2 5') 
cct.add('C1 2 0 1') 


pprint(cct.V)

pprint(cct.I)




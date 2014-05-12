from lcapy import pprint, Circuit

cct = Circuit('Transformer')

cct.net_add('Vs 1 0 10') 
# Cannot have floating voltage source so tie to ground with arbitrary R
cct.net_add('R1 2 0 100') 
cct.net_add('TF1 2 0 1 0 10') 



pprint(cct.V)

pprint(cct.I)




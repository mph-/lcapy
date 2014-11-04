from lcapy import pprint, Circuit

cct = Circuit('Transformer')

cct.add('Vs 1 0 dc 10') 
# Cannot have floating voltage source so tie to ground with arbitrary R
cct.add('R1 2 0 100') 
cct.add('TF1 2 0 1 0 10') 



pprint(cct.V)

pprint(cct.I)




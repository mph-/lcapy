from lcapy import pprint, Circuit

cct = Circuit()
cct.add('V1 1 0 step')
cct.add('R1 1 2')
cct.add('L1 2 0')
cct.add('L2 3 0')
cct.add('K1 L1 L2 1')
cct.add('R2 3 0')


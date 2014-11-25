from lcapy import pprint, Circuit

cct = Circuit()

cct.add('V_s fred 0 step') 
cct.add('R_a fred bert') 
cct.add('R_b bert 0') 


pprint(cct.V)

pprint(cct.I)




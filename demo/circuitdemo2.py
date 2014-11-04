from lcapy import pprint, Circuit

cct = Circuit('Voltage divider')

cct.add('V_s fred 0 dc') 
cct.add('R_a fred bert') 
cct.add('R_b bert 0') 


pprint(cct.V)

pprint(cct.I)




from lcapy import pprint, Circuit

cct = Circuit('Voltage divider')

cct.net_add('V_s fred 0 dc') 
cct.net_add('R_a fred bert') 
cct.net_add('R_b bert 0') 


pprint(cct.V)

pprint(cct.I)




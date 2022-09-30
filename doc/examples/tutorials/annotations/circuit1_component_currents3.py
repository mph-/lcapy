from lcapy import *

cct = Circuit('circuit1.sch')
cct.annotate_currents(('R1', 'R2', 'R3', 'R4'), evalf=False, flow=True).draw(__file__.replace('.py', '.png'), draw_nodes='connections')



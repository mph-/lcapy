from lcapy import *

cct = Circuit('circuit1.sch')
cct.annotate_voltages(('R1', 'R2', 'R3', 'R4')).draw(__file__.replace('.py', '.png'), draw_nodes='connections')



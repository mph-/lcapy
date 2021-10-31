from lcapy import *

cct = Circuit('circuit1.sch')
cct.annotate_node_voltages().draw(__file__.replace('.py', '.png'), draw_nodes='primary')



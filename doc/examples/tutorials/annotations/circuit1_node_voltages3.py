from lcapy import *

cct = Circuit('circuit1.sch')
cct.annotate_node_voltages(label_voltages=True, show_units=False).draw(__file__.replace('.py', '.png'), draw_nodes='primary', node_spacing=2.5)



from lcapy import Circuit
from lcapy.circuitgraph import CircuitGraph

a = Circuit('graph2.sch')
G = CircuitGraph(a)
G.draw(__file__.replace('.py', '.png'))

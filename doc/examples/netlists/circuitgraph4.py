from lcapy import Circuit
from lcapy.circuitgraph import CircuitGraph

a = Circuit('graph4.sch')
G = CircuitGraph(a)
G.draw(__file__.replace('.py', '.png'))

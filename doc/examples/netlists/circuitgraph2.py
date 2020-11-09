from lcapy import *

a = Circuit('graph2.sch')
G = CircuitGraph(a)
G.draw(__file__.replace('.py', '.png'))

   

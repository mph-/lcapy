from lcapy import *

a = Circuit('graph4.sch')
G = CircuitGraph(a)
G.draw(__file__.replace('.py', '.png'))

   

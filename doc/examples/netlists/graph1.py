from lcapy import *

a = Circuit("""
V1 1 0 {1 + u(t)}
R1 1 2
L1 2 0
R2 1 3
L2 3 0""")
G = CircuitGraph(a)
G.draw(__file__.replace('.py', '.png'))

   

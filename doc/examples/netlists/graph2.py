from lcapy import *

a = Circuit("""
V1 1 0; down
R1 1 2; right
L1 2 3; right
R2 3 4; right
L2 2 0_2; down
C2 3 0_3; down
R3 4 0_4; down
W 0 0_2; right
W 0_2 0_3; right
W 0_3 0_4; right""")
G = CircuitGraph(a)
G.draw(__file__.replace('.py', '.png'))

   

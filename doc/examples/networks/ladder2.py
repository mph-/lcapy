from lcapy import R, C
n = C(2) | (R(3) + (C(4) | (R(5) + (C(6)))))
n.draw(__file__.replace('.py', '.png'), layout='ladder')


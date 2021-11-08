from lcapy import R, L
n = L(2) | (R(3) + (L(4) | (R(5) + (L(6)))))
n.draw(__file__.replace('.py', '.png'), layout='ladder')


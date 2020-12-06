from lcapy import R, C, ladder
n = ladder(C(1), R(2), C(3), R(4))
n.draw(__file__.replace('.py', '.png'), form='ladder')


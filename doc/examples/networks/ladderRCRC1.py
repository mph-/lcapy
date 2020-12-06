from lcapy import R, C, ladder
n = ladder(R(1), C(2), R(3), C(4))
n.draw(__file__.replace('.py', '.png'), form='ladder')


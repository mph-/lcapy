from lcapy import R, C
n = R('R1') + (C('C2') | (R('R2') + (C('C3') | R('R3'))))
n.draw(__file__.replace('.py', '.png'), layout='ladder')


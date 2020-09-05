from lcapy import R, C
n = C('C1') | (R('R1') + (C('C2') | (R('R2') + (C('C3') | (R('R3') + C('C4'))))))
n.draw(__file__.replace('.py', '.png'), form='ladder')


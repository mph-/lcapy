from lcapy import R, C, L
n = C('C1') | ((R('R1') + L('L1')) + (C('C2') | ((R('R2') + L('L2')) + (C('C3') | (R('R3') + L('L3')) + C('C4')))))
n.draw(__file__.replace('.py', '.png'), layout='ladder')


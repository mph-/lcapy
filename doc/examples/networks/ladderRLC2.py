from lcapy import R, C, L
n = ((R('R1') + L('L1')) + (C('C2') | ((R('R2') + L('L2')) + (C('C3') | (R('R3') + L('L3'))))))
n.draw(__file__.replace('.py', '.png'), layout='ladder')


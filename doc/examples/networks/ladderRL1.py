from lcapy import R, L
n = L('L1') | (R('R1') + (L('L2') | (R('R2') + (L('L3') | R('R3')))))
n.draw(__file__.replace('.py', '.png'), layout='ladder')


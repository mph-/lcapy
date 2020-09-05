from lcapy import R, L
n = (R('R1') | L('L1')) + (R('R2') | L('L2'))
n.draw(__file__.replace('.py', '.png'))


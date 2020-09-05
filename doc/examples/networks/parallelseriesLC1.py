from lcapy import C, L
n = (C('C1') + L('L1')) | (C('C2') + L('L2'))
n.draw(__file__.replace('.py', '.png'))


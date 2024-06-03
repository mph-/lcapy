from lcapy import R, L, C

n = R('R0') + (C('C1') | (R('R1') + L('L1')))
n.draw('RRLC.png')

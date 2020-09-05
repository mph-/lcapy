from lcapy import R, C
n = (R('R1') | C('C1')) + (R('R2') | C('C2'))
n.draw(__file__.replace('.py', '.png'))


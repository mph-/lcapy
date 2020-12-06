from lcapy import R, C, ladder
n = C('C1') | (R('R1') + (C('C2') | (R('R2') + (C('C3')))))
n2 = ladder(None, C('C1'), R('R1'), C('C2'), R('R2'), C('C3'))
n.draw(__file__.replace('.py', '.png'), form='ladder')


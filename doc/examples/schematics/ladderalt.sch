# print(Ladder(Z('Z'),Y('Y'),Z('Z'),Y('Y')).netlist())
P1 1 0; down, v=V_1
W 1 9; right=0.5, i>_=I_1
W 0 8; right
Y1 9 8 Y; down, l=OP1
W 9 7; right=0.5
W 8 6; right=0.5
Z2 7 5 Z; right, l=OP2
W 6 4; right
O 7 6; down
O 5 4; down
Y2 13 12 Y; down, l=OP3
W 5 13; right=0.5
W 13 3; right=0.75, i<_=I_2
W 4 12; right=0.5
W 12 2; right=0.5
P 3 2; down, v=V_2
;label_nodes=False, draw_nodes=connections

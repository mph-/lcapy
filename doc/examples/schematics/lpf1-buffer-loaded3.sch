Vi 1 0_1; down
Rs 1 2; right=1.5
C 2 0; down
W 0_1 0; right
W 0 0_2; right
Rin 2_2 0_2; down, v=V_{in}
W 2 2_2; right
E1 3 0_3 2 0 A; down, l=A V_{in}
Rout 3 4; right=1.5
RL 4 0_4; down, v=V_o
W 0_2 0_3; size=1.2
W 0_3 0_4
P1 4 0_4; down
; style=european
C1 1 0_1; down, dashed, blue
V1 2 1; down
W 2 3; up
W 3 4; right=1.5
V2 5 4; down=1.25, fixed, l=V_d/2
V3 4 22; down=1.25, fixed, l=V_d/2
W 5 6; right
W 6 7; down
W 7 8; right=2
W 22 21; right
W 21 20; up
W 20 19; right=2
TF 8 9 19 18 core; up, l={}
W 9 10; right=2
W 18 17; right=2
W 10 11; up
W 11 12; right
W 17 16; down
W 16 15; right
R1 12 13; down
R2 13 15; down
# Hack
O 4 13; right
W 13 14; right=1.5
W 14 23; down
C2 23 0_23; down, dashed, blue
W 0_1 0_23; right
W 0 0_1; right=0.25
W 0_23 0_24; right=0.25
; draw_nodes=connections


V1 1 0_1; down
R1 1 2 R_s; right=1.5
C1 2 0; down, v=v_C
W 0_1 0; right
W 0 0_2; right
R2 2_2 0_2 R_in; down, v=v_{in}
W 2 2_2; right
E1 3 0_3 2 0 A; down, l=A v_{in}
R3 3 4 R_out; right=1.5
R4 4 0_4 R_L; down
W 0_2 0_3; size=1.2
W 0_3 0_4
P1 4 0_4; down
; label_ids=False, draw_nodes=connections, label_nodes=False




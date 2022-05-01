R1 1 2; right
L1 2 3 {L_1 - a * M}; right
L3 4 0_4 {a * M}; down
P1 1 0; down, v_=v_1
W 0 0_3; right
W 3 4; right=0.5
W 0_3 0_4; right=0.5
W 4 5; right=0.5
W 0_4 0_7; right=0.5
L2 5 6 {a**2 * L_2 - a * M}; right
R2 6 7 {a**2 * R_2}; right
W 7 10; right=0.5
TF 8 0_8 10 0_7; right, l={N_1:N_2}
W 8 9; right
W 0_8 0_9; right
P2 9 0_9; down, v^=v_2
;label_nodes=False, draw_nodes=connections, label_ids=False

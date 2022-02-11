R1 1 2; right
L1 2 3 {L_1 - a * M}; right
L3 3 0_3 {a * M}; down=1.5
P1 1 0; down, v=v_1
W 0 0_3; right
W 3 4; right
W 0_3 0_4; right
TF 5 0_5 4 0_4; right, l={N_1:N_2}
W 5 6; right=0.5
W 0_5 0_6; right=0.5
L2 6 7 {L_2 - M / a}; right
R2 7 8; right
P2 8 0_8; down, v^=v_2
W 0_6 0_8; right
;label_nodes=False, draw_nodes=connections, label_ids=False

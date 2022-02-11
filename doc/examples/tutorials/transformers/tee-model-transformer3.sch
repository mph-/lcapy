L1 1 3 {L_1 - M * a}; right, i>^=i_1
L3 3 0_3 {M * a}; down=1.5
P1 1 0; down, v_=v_1
W 0 0_3; right
W 0_3 0_4; right
L2 3 6 {L_2 *a**2 - M * a}; right
W 6 4; right=0.5
TF 5 0_5 4 0_4 {1/a}; right
W 5 7; right=0.5, i^<=i_2
W 0_5 0_7; right=0.5
P2 7 0_7; down, v^=v_2
;label_nodes=False, draw_nodes=connections, label_ids=False

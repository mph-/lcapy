Vs 6 0_3 ac; down=2
Cs 6 4; right
W 0_3 0; right
E1 1 0 opamp 0_4 2_1 A; right, mirror
W 1 1_2; right=0.5
P 1_2 0_6; down
W 4 4_3; right
Vn 4_3 2 noise; right
W 2 2_1; right=0.5
In 2 0_1 noise; down
W 0_4 0_2; down
W 0 0_5; right
W 0_5 0_1; right
W 0_1 0_2; right
W 0_2 0_6; right=0.5
W 4 4_2; up
VnR 4_2 3 noise {sqrt(4 * k_b * T * R)}; right
R 3 1_1; right
W 1_1 1; down
; draw_nodes=connections, label_nodes=primary

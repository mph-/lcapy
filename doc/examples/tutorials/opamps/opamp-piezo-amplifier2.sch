Cs 1 0; down=4
W 1 1_1; right
Rs 1_1 0_1; down=4
W 0 0_1; right
W 1_1 1_2; right=2
Vn 1_2 2 noise; right
E 5_1 0 opamp 2 3_2 A; right
W 3 3_1; right
W 3_1 3_2; right
W 3 3_3; down
R1 3_3 4; down
W 4 0_2; down
W 0_1 0_2; right
W 3_1 3_4; down=1.5
Inn 3_4 0_3 noise; down
W 0_2 0_3; right
W 2 2_1; down=2
Inp 2_1 0_4 noise; down
W 0_3 0_4; right
W 3_3 3_5; right=3
R2 3_5 5_2; right
W 5_1 5_2; down=2
W 5_1 5; right
Po 5 0_5; down, v=v_o
W 0_4 0_5; right
; draw_nodes=connections, label_nodes=none

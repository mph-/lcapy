Vsp 1 0_1 ; down
W 0_1 0; down=0.1, 0V
W 1 1_1; right
R1 1_1 6; right
Vsm 2 0_3 ; down
W 0_3 0; down=0.05, 0V
R3 2 3; right
E1 5_2 4_2 fdopamp 3 6 0_4 A; right, mirror, l
W 5_2 5; right
W 4_2 4; right
P2 5 4; down
W 5_1 5_2; down
W 6_1 6; down
W 3 3_1; down
R2 6_1 5_1; right
W 4_2 4_1; down
R4 3_1 4_1; right
W 0 0_4; right=0.4
W 0 0_2; down=0.1, 0V
; draw_nodes=connected

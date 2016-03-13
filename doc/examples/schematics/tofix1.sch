R1 1 b {Rs / 2}; right
W b 2; right=0.5
W 2 2_1; up=0.5
W 2 2_2; down=0.5
R2 2_1 3_1 Rp; right
C1 2_2 3_2 Cp; right
W 3 3_1; up=0.5
W 3 3_2; down=0.5
W 3 4; right=0.5
W 4 4_1; up=0.5
W 4 4_2; down=0.5
C2 4_2 5_2 Cw; right
R3 4_1 5_1 Rw; right
W 5 5_1; up=0.5
W 5 5_2; down=0.5
W 5 6; right=0.5
W 6 6_1; up=0.5
W 6 6_2; down=0.5
C3 6_2 7_2 Cp; right
R4 6_1 7_1 Rp; right
W 7 7_1; up=0.5
W 7 7_2; down=0.5
W 7 c; right=0.5
R5 c 8 {Rs / 2}; right
W b 2_3; up
W c 7_3; up=1.25
R6 2_3 7_3 Rb; right
# FIXME, the following open circuit should not be required
#O 2 3; right
; draw_nodes=connections, label_nodes=false, label_ids=false


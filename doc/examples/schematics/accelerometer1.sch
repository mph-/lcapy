Pm1 1 0_1; down, v=f_1
W 1 0_1; down, dashed
Pm2 2 0_2; down, v=f_2
W 2 0_2; down, dashed
Lm1 1 3; right=3, i>^=u_1
Lm2 3 2; right=2, i^<=u_2
Rm 3 4; down
Cm 4 5; down
TF 5 6 9 0 k; right
W 6 0_6; down
W 0_1 0_10; right
W 0_10 0_6; right, free
W 0_6 0_2; right

Pe 7 0_7; down, v=v
W 7 8; right=0.5, i=i
W 0_7 0_8; right=0.5
W 8 9; right
W 0_8 0; right
C0 8 0_8; down
W 0 0_10; down, dashed
; label_nodes=none, draw_nodes=connections

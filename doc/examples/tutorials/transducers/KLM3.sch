Pb 1 0_1; down, v_=F_b
TLb 2 6 1 0_1 Z_0 {s / c} {d / 2}; right

Pf 3 0_3; down, v^=F_f
TLf 3 0_3 15 5 Z_0 {s / c} {d / 2}; right

W 2 15; right=0.5
W 2 4; down=1.5
W 6 5; right=0.5

TF1 4 7 8 9 phi; right
W 9 7; right=0.5, dashed
W 7 10; right=0.5
W 5 10; down

Pe 11 0; down, v_=V
C0 11 12; right, i>^=I
Z 12 8 {j * X_1}; right
W 0 9; right

W 9 7; right, ignore

; draw_nodes=connections, label_nodes=none, label_ids=none

Rs 1 _anon1; down
VnRs _anon1 0 noise {sqrt(4 * k * T * Rs)}; down
Vn 1 2 noise; right
W 2 3; right
In1 2 0_2 noise; down, l=I_{n+}
W 0 0_2; right
In2 5 0_5 noise; down, l=I_{n-}
W 5 4; right
W 0_2 0_5; right
W 4 6; down
R1 6 _anon2; down
VnR1 _anon2 0_6 noise {sqrt(4 * k * T * R1)}; down
W 0_5 0_6; right
R2 6 _anon3; right
VnR2 _anon3 7 noise {sqrt(4 * k * T * R2)}; right
W 8 7; down
E 8 0 opamp 3 4 A; right
W 8 9; right
W 0_6 0_9; right
P 9 0_9; down
; draw_nodes=connections, label_nodes=none
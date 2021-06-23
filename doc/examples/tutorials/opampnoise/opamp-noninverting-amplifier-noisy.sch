NRs 1 _nodeanon1 Rs; down
VnRs _nodeanon1 0 noise {sqrt(4 * k_B * T * Rs)}; down
Vn 1 2 noise; right
W 2 3; right
In1 2 0_2 noise; down, l=I_{n+}
W 0 0_2; right
In2 5 0_5 noise; down, l=I_{n-}
W 5 4; right
W 0_2 0_5; right
W 4 6; down
NR1 6 _nodeanon2 R1; down
VnR1 _nodeanon2 0_6 noise {sqrt(4 * k_B * T * R1)}; down
W 0_5 0_6; right
NR2 6 _nodeanon3 R2; right
VnR2 _nodeanon3 7 noise {sqrt(4 * k_B * T * R2)}; right
W 8 7; down
E 8 0 opamp 3 4 A; right
W 8 9; right
W 0_6 0_9; right
P 9 0_9; down
; draw_nodes=connections, label_nodes=none
P PIO GND; down
W PIO 1; right
W GND 3; right
D1 1 2; up
D2 3 1; up
W 1 4; right
W 4 5; up, size=0.25
W 4 6; down, size=0.25
M1 7 5 8 pmos; right
M2 7 6 9 nmos; right
W 3 9; right
W 2 8; right
W 7 10; right, size=0.5
W 8 11; right
W 9 12; right
VDD 11 12;down
; draw_nodes=connections, label_nodes=alpha

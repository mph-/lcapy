P 1 2; down, v=V_1
W 1 5; right=0.5, i=I_1
W 2 6; right=0.5, i<=I_1
TPB 7 8 5 6; right, l={Source-free B-parameters two-port}
V2b 9 7; left
W 8 10; right
I2b 10 9; up
W 9 3; right=0.5, i=I_2
W 10 4; right=0.5, i<=I_2
P 3 4; down, v^=V_2
; label_nodes=none, draw_nodes=connections

P 1 2; down, v_=V_1
W 1 11; right=0.5, i=I_1
W 2 12; right=0.5, i<=I_1
I1g 12 11; up
W 11 5; right=0.25
W 12 6; right=0.25
TPG 7 8 5 6; right, l={Source-free G-parameters two-port}
V2g 9 7; left
W 8 10; right
W 9 3; right=0.5, i=I_2
W 10 4; right=0.5, i<=I_2
P 3 4; down, v^=V_2
; label_nodes=none, draw_nodes=connections

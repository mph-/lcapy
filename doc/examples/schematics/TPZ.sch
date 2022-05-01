P 1 2; down, v_=V_1
W 1 11; right=0.25, i=I_1
W 2 12; right=0.25, i<=I_1
V1z 11 5; right
W 12 6; right
TPZ 7 8 5 6; right, l={Source-free Z-parameters two-port}
V2z 9 7; left
W 8 10; right
W 9 3; right=0.25, i=I_2
W 10 4; right=0.25, i<=I_2
P 3 4; down, v^=V_2
; label_nodes=none, draw_nodes=connections

W 1 VDD; up=0.05, implicit, l=V_{DD}
Rd 1 3a; down, f={I_{ds} + I_{o}}
W 3a 3; down=0.5
M1 3 4 5; right, kind=nfet
W 3 o; right, f=I_o
P o 7; down, v=V_o
W 7 0; down=0.05, implicit, l=0V
Rs 5 5a; down, f=I_{ds}
W 5a VSS; down=0.05, implicit, l=0V
Rl o 9; right
W 9 6; right=0.5
Vl 6 8; down
W 8 0; down=0.05, implicit, l=0V
W 11 4; right=0.5
P 11 13; down, v=V_{i}
W 13 0; down=0.05, implicit, l=0V
; draw_nodes=connections, label_nodes=none, bipole voltage style={color=blue}, bipole flow style={color=blue}

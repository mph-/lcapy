W 1a VDD; up=0.05, implicit, l=V_{DD}
W 1a 1; down=0.5, f={I_{dn}}
M1 1 2 3a; right, kind=nfet
M2 5 4 3b; right, kind=pfet
W 3a 3; down=0.5
W 3 3b; down=0.5
W 3 o; right, f=I_o
P o 7; down, v=V_o
W 7 0; down=0.05, implicit, l=0V
W 5 5a; down=0.5, f<=I_{dp}
W 5a VSS; down=0.05, implicit, l=V_{SS}
Rl o 9; right
W 9 6; right=0.5
Vl 6 8; down
W 8 0; down=0.05, implicit, l=0V
W 10 2; right=0.5
W 11 4; right=0.5
P 10 12; down, v=V_{gn}, bipole voltage style={color=blue}
P 11 13; down, v=V_{gp}
W 12 0; down=0.05, implicit, l=0V
W 13 0; down=0.05, implicit, l=0V

; draw_nodes=connections, label_nodes=none, bipole voltage style={color=blue}, bipole flow style={color=blue}

U1 inverter ._in ._VSS ._out ._VDD; right
W 5 U1._in; right=0.5
W U1._out 6; right=0.5
W 6 9; right=0.5
R1 7 8; right
W 5 1; down=0.75
W 6 2; down=0.75
W 5 7; up=0.75
W 6 8; up=0.75
XT1 1 2; right
C1  1 4; down=0.8
C2  2 3; down=0.8
W 4 0; down=0.1, implicit, l=GND
W 3 0; down=0.1, implicit, l=GND
; draw_nodes=connections, label_nodes=False

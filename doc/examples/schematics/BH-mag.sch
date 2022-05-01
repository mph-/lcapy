VHV 2c 1c; down, color=blue
D1 2c 3; right, color=blue
W 3 3a; right=0.5, color=blue
W 1c 1b; right, color=blue
W 1b 1d; right=0.5, color=blue
C 3a 1d; down, color=blue
W 4 3; down=0.5, color=blue
W 4 4a; right=2, color=blue
MC1 4a 5 6; right=1.5, scale=1.5, bodydiode, kind=nfet, l=QC1, color=blue
W 1b 1; down=0.5, color=blue
W 1 1a; right
MC2 6 7 1a; right=1.5, scale=1.5, bodydiode, kind=nfet, l=QC2
MD1 1a 8 9; right=1.5, scale=1.5, bodydiode, kind=nfet, l=QD1
MD2 9 16 0a; right=1.5, scale=1.5, bodydiode, kind=nfet, l=, color=blue
W 1 1e; down=0.5, color=blue
VLV 1e 0; down, color=blue
W 0 0a; right, color=blue
W 6 6a; right=2, color=blue
W 6a 6b; right=2
MA1 6a 10 12; right=1.5, scale=1.5, bodydiode, kind=nfet, l=QA1, color=blue
MA2 12 11 9a; right=1.5, scale=1.5, bodydiode, kind=nfet, l=QA2
W 9 9a; right=2, color=blue
MB1 6b 14 13a; right=1.5, scale=1.5, bodydiode, kind=nfet, l=QB1
MB2 13a 15 9b; right=1.5, scale=1.5, bodydiode, kind=nfet, l=QB2, color=blue
W 9a 9b; right, color=blue
W 12 12a; right=0.5, color=blue
L 12a 13; right, i>_=i_L, v=v_L, color=blue
W 13 13a; right=0.5, color=blue
# Hack for Circuitikz bounding box problem
O 13a 13b; right=0.8
; draw_nodes=connections, label_nodes=none

C1 1 0_1; down, dashed, blue
V1 2 1; down, l=v_c
W 2 4; up=0.5, i=i_c
V2 5 4; down=1.25, fixed, l=v_d
W 5 6; right
W 6 7; down
W 7 8; right=2, i={(1 - a) i_c + i_d}
W 4 3; right=1, fixed
W 3 19; right, i={a i_c - i_d}
TF 8 9 19 18 core; up, l={}
W 9 10; right=2
W 18 13; right=2
W 10 11; up
W 11 12; right
R1 12 13; down, l=R_L
# Hack
O 4 13; right
O 6 11; right
W 13 23; down, i=i_c
C2 23 0_23; down, dashed, blue
W 0_1 0_23; right, l={common ground}
W 0 0_1; right=0.25
W 0_23 0_24; right=0.25
;;\node[blue,draw,dashed,inner sep=5mm,anchor=north, fit=(1) (5) (7), label=Device 1] {};
;;\node[blue,draw,dashed,inner sep=5mm,anchor=north, fit=(11) (13), label=Device 2] {};
;draw_nodes=connections, label_nodes=false



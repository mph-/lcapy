U1 inverter ._IN ._VSS ._OUT ._VDD; right, size=2, scale=2, l={}
U2 inverter ._IN ._VSS ._OUT ._VDD; right, size=2, scale=2, l={}
W U1._OUT U2._IN; right=1
W U1._VSS 4_2; down=0.3, implicit, l=0V
W U2._VSS 4_3; down=0.3, implicit, l=0V
W U1._VDD 3_2; up=0.3, implicit, l=V_{DD1}
W U2._VDD 3_3; up=0.3, implicit, l=V_{DD2}
W U2._IN 7_1; right=0.25, dashed, fixed, color=red
D 7_1 6_1; up=0.25, scale=0.5, l=, color=red
D 6_2 7_1; up=0.25, scale=0.5, l=
W 6_1 U2._VDD; right=0.6, dashed, color=red
W 6_2 U2._VSS; right=0.6, dashed
M1 1 2 3 pmos; right=0.25, scale=0.25, l={}, color=red
M2 1 4 5 nmos; right=0.25, scale=0.25, l={}
W 3 3_1; up=0.0, color=red
W 3_1 U1._VDD; right=0.2, dashed, color=red
W 5 5_1; down=0.0
W 5_1 U1._VSS; right=0.2, dashed
W 1 U1._OUT; right, dashed, color=red
;draw_nodes=connections,label_nodes=false,thickness=2


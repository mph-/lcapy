U1 inverter; right=2, l={}
U2 inverter; right=2, l={}
W U1.out U2.in; right=1, color=red
W U1.vss 4_2; down=0.3, implicit, l=0V
W U2.vss 4_3; down=0.3, implicit, l=0V
W U1.vdd 3_2; up=0.3, implicit, l=V_{DD1}, color=red
W U2.vdd 3_3; up=0.3, implicit, l=V_{DD2}, color=red
W U2.in 7_1; right=0.25, dashed, fixed, color=red
D 7_1 6_1; up=0.25, scale=0.5, l=, color=red
D 6_2 7_1; up=0.25, scale=0.5, l=
W 6_1 U2.vdd; right=0.6, dashed, color=red
W 6_2 U2.vss; right=0.6, dashed
M1 1 2 3 pmos; right=0.25, scale=0.25, l={}, color=red
M2 1 4 5 nmos; right=0.25, scale=0.25, l={}
W 3 3_1; up=0.0, color=red
W 3_1 U1.vdd; right=0.2, dashed, color=red
W 5 5_1; down=0.0
W 5_1 U1.vss; right=0.2, dashed
W 1 U1.out; right, dashed, color=red
;draw_nodes=connections,label_nodes=false,thickness=2


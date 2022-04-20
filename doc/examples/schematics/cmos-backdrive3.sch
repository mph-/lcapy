# This has a stretch conflict, it needs a wire.
U1 inverter; right=2, l={}
W U1.vss 4_2; down=0.3, implicit, l=0V
W U1.vdd 3_2; up=0.3, implicit, l=V_{DD1}
M1 1 2 3 pmos; right=0.25, scale=0.25, l={}
M2 1 4 5 nmos; right=0.25, scale=0.25, l={}
W 3 U1.vdd; right=0.4, dashed
W 5 U1.vss; right=0.4, dashed
W 1 U1.out; right, dashed
;draw_nodes=connections,label_nodes=false,thickness=2

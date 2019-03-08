U1 buffer; right, pinnames={vdd=VDDIO}
W U1.out 1; right=0.5
R 1 2; right=2, i>^=I_o
D 2 0_2 led; down, v=V_f, l={}
W U1.vss 0_3; down
W 0_3 0_1; right
W 0_1 0_2; right
O 1 0_1; down, v=V_o
W 0_3 0; down=0.1, implicit, l=0V
W U1.vdd _vddIO; up=0.3, implicit, l=V_{DDIO}
; draw_nodes=connections, label_nodes=false

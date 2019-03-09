U1 chip4141; right, pinlabels={out1=PIO1,out3=PIO2}
W U1.vdd VDD; implicit, up=0.2, l=3V3
W U1.vss 0; implicit, down=0.5, l=0V
R1 U1.out3 1; right
D1 1 3 led; down
W 3 0; down=0.2, implicit, l=0V
R2 U1.out1 2; right
W 2 5; right
D2 5 4 led; down
W 4 0; down=0.7, implicit, l=0V
; draw_nodes=connections, help_lines=1

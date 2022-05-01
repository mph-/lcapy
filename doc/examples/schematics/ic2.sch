U1 chip4141; right, pinlabels={r1=PIO1,r3=PIO2}
W U1.vdd VDD; implicit, up=0.2, l=3V3
W U1.vss 0; implicit, 0V
R1 U1.r3 1; right
D1 1 3 led; down
W 3 0; down=0.2, 0V
R2 U1.r1 2; right
W 2 5; right
D2 5 4 led; down
W 4 0; down=0.7, 0V
; draw_nodes=connections, help_lines=1

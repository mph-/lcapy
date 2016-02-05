#U1 chip4141 _1 _2 _3 _4 Vss _5 PIO1 PIO2 _8 Vdd; right
U1 chip4141 1 2 3 4 Vss 5 PIO1 PIO2 8 Vdd; right
W U1.Vdd Vdd; implicit, up=0.2, l=3V3
W U1.Vss 0; implicit, down=0.5, l=0V
R1 U1.PIO1 1; right
D1 1 3 led; down
W 3 0; down=0.2, implicit, l=0V
R2 U1.PIO2 2; right
W 2 5; right
D2 5 4 led; down
W 4 0; down=0.7, implicit, l=0V
; draw_nodes=all, help_lines=1




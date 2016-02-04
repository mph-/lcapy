U1 chip4141 _1 _2 _3 _4 Vss _5 PIO1 PIO2 _8 Vdd; right
W U1.Vdd Vdd; implicit, up, l=3V3
W U1.Vss 0; implicit, down, l=0V
R1 U1.PIO1 1; right
D1 1 3 led; down
W 3 0; down, implicit, l=0V
R2 U1.PIO2 2; right
W 2 5; right
D2 5 4 led; down
W 4 0; down, implicit, l=0V, size=1
; draw_nodes=connections


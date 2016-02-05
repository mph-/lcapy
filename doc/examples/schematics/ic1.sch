; draw_nodes=all, help_lines=1
U1 chip2121 _1 _2 Vss PIO1 PIO2 Vdd; right
W U1.Vdd Vdd; implicit, up=0.2, l=3V3
W U1.Vss 0; implicit, down=0.8, l=0V
R1 U1.PIO1 1; right
D1 1 3 led; down
W 3 0; down=0.1, implicit, l=0V
R2 U1.PIO2 2; right
W 2 5; right
D2 5 4 led; down
W 4 0; down=0.1, implicit, l=0V, size=1



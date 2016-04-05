; draw_nodes=connections, help_lines=1
U1 chip2121 _1 _2 .VSS .PIO1 .PIO2 .VDD; right=2, pins=auto
W U1.VDD VDD; implicit, up=0.2, l=3V3
W U1.VSS 0; implicit, down=0.7, l=0V
R1 U1.PIO1 1; right
D1 1 3 led; down
W 3 0; down=0.1, implicit, l=0V
R2 U1.PIO2 2; right
W 2 5; right
D2 5 4 led; down
W 4 0; down=0.1, implicit, l=0V




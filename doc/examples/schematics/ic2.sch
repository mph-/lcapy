U1 chip4141 .1 .2 .3 .4 .VSS .5 .PIO1 .PIO2 .8 .VDD; right, pins={,,,,VSS,,PIO1,PIO2,,VDD}
W U1.VDD VDD; implicit, up=0.2, l=3V3
W U1.VSS 0; implicit, down=0.5, l=0V
R1 U1.PIO1 1; right
D1 1 3 led; down
W 3 0; down=0.2, implicit, l=0V
R2 U1.PIO2 2; right
W 2 5; right
D2 5 4 led; down
W 4 0; down=0.7, implicit, l=0V
; draw_nodes=connections, help_lines=1
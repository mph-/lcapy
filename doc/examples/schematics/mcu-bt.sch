U1 chip3131 ._UNUSED1 ._UNUSED2 ._UNUSED3 .VSS .RXD .TXD .PIO1 .VDD; right=2, l=MCU, pins=auto, aspect=0.75
U2 chip3131 ._UNUSED1 .RXD .TXD .VSS ._UNUSED2 ._UNUSED3 ._UNUSED4 .VDD; right=2, l=Bluetooth, pins=auto, aspect=0.75
W U1.TXD U2.RXD; right=2.5
W U1.RXD U2.TXD; right=2.5
U3 chip1310 ._IN .EN .GND ._NC ._OUT; right=1.5, aspect=1.5, pins=auto
W U1.PIO1 1; right=0.25
W 1 U3.EN; up
W U3._OUT 2; right=0.5
W 2 U2.VDD; down
U4 chip1310 ._IN ._EN .GND ._NC ._OUT; right=1.5, aspect=1.5, pins=auto
W U4._OUT 3; right=0.5
W 3 U1.VDD; down
# Hack
#O 3 U3._IN; right
W 5 U4._IN; right=0.5
W 5 _5; up=0.4, implicit, l=V_{bat}
W 6 U3._IN; right=0.5
W 6 _6; up=0.4, implicit, l=V_{bat}
W U1.VSS 0_5; down=0.1, implicit, l=0V
W U2.VSS 0_6; down=0.1, implicit, l=0V
W U4.GND 0_7; down=0.1, implicit, l=0V
W U3.GND 0_8; down=0.1, implicit, l=0V
R 1 7; right
W 7 0_9; down=0.1, implicit, l=0V
; draw_nodes=connections, label_nodes=alpha, node_spacing=2
U1 chip3131 ._UNUSED1 ._UNUSED2 ._UNUSED3 .VSS .RXD .TXD .PIO .VDD; right, l=MCU
U2 chip2121 .RXD .TXD .VSS ._UNUSED1 ._UNUSED2 .VDD; right, l=Bluetooth
W U1.TXD U2.RXD; right=2
W U1.RXD U2.TXD; right=2
U3 chip1310 ._IN .EN .GND ._NC ._OUT; right
W U1.PIO 1; right=0.5
W 1 U3.EN; up
W U3._OUT 2; right
W 2 U2.VDD; down
U4 chip1310 ._IN .EN .GND ._NC ._OUT; right
W U4._OUT 3; right
W 3 U1.VDD; down
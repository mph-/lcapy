U1 chip2121 ._UNUSED1 ._UNUSED2 .VSS .PWM2 .PWM1 .VDD; right, l=MCU, pins={,,VSS, PWM2, PWM1, VDD}
M1 9 10 11; right
M2 12 13 14; right
W U1.PWM1 10; right=0.1
W U1.PWM2 13; right=0.1
TF1 1 2 3 4 tapcore 5 _6; right, l=1:100
W 9 8; up=0.1
W 8 3; right
W 7 5; right=0.5
W 12 4; up=0.1
W 14 14_1; down=0.1, implicit, l=0V
W 11 11_1; down=1.1, implicit, l=0V
W 7 7_1; up=0.8, implicit, l=3.9V
W U1.VSS _16; down=0.4, implicit, l=0V
W U1.VDD _17; up=0.2, implicit, l=3.3V
W 1 1_1; right=0.5
W 2 2_1; right=0.5
; draw_nodes=connections, label_nodes=alpha






U1 chip3131; right=1.8, l=MCU, pins={out3=RXD, out2=TXD, out1=PIO1}, aspect=0.8
U2 chip3131; right=1.8, l={\hspace{5mm}Bluetooth}, pins={in2=RXD, in3=TXD}, aspect=0.8
W U1.out2 U2.in2; right=2.5
W U1.out3 U2.in3; right=2.5
U3 regulator; right=1.5, aspect=1.5, pins={en=EN}, l=VREG
W U1.out1 1; right=0.25
W 1 U3.en; up
W U3.out 2; right=0.5
W 2 U2.vdd; down
U4 regulator; right=1.5, aspect=1.5, l=VREG
W U4.out 3; right=0.5
W 3 U1.vdd; down
# Hack
O 3 U3.in; right
W 5 U4.in; right=0.5
W 5 _5; up=0.4, implicit, l=V_{bat}
W 6 U3.in; right=0.5
W 6 _6; up=0.4, implicit, l=V_{bat}
W U1.vss 0_5; down=0.1, implicit, l=0V
W U2.vss 0_6; down=0.1, implicit, l=0V
W U4.gnd 0_7; down=0.1, implicit, l=0V
W U3.gnd 0_8; down=0.1, implicit, l=0V
R 1 7; right
W 7 0_9; down=0.1, implicit, l=0V
; draw_nodes=connections, label_nodes=pins, node_spacing=2

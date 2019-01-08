W PIO _2;right, size=0.5
W GND _9;right, size=0.5
P PIO GND;down
W _2 _3;up, size=0.5
W _2 _4;down, size=0.5
W _5 _6;up, size=0.5
W _5 _7;down, size=0.5
R _3 _6 R;right
C _4 _7 C;right
W _5 _8;right
VDD _8 _9 {0.5 * V_DD};down
;draw_nodes=connections, label_ids=False
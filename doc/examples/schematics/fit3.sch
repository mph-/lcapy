U1 chip3333; right, l={MCU}, pinlabels={vss,vdd}
R1 1 1a; down
V_DD 1a 0; down
L1 1 2; right=1
L2 0 0_2; right
W 2 5; right=0.5
W 0_2 0_5; right=0.5
Cbulk 5 0_5; down=2, kind=electrolytic
L3 5 6; right
L4 0_5 0_6; right
W 6 3; right=0.5
W 0_6 0_3; right=0.5
L5 3 4; right, color=blue
L6 0_3 0_4; right, color=blue
W 4 U1.vdd; down=0.25, color=blue
W 0_4 U1.vss; up=0.25, color=blue
Clocal 6 0_6; down, color=blue
A U1.tl
A U1.br
;;\node[blue,draw,dashed,inner sep=8mm, fit=(Cbulk) (U1@tl) (U1@br), label=PCB]{};
; label_nodes=none, draw_nodes=connections, label_values=false, label_ids=false

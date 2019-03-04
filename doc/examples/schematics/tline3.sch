U1 buffer; right
TL1 1 2 U1.out 4; right=2, scale=2, l={50\,$\Omega$}
R 1 2; down
W U1.vdd 9; up=0.1, implicit, l=Vdd
W U1.vss 11; down
W 11 4; right=0.5
W 11 10; down=0.1, implicit, l=Vss
W 5 U1.in; right=0.5
;draw_nodes=connections, help_lines=0.1

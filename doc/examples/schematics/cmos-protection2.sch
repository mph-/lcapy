U1 inverter ._IN ._VSS ._OUT ._VDD; right
R 2 U1._IN; right
D1 2 3; up
W 3 3_1; up=0.1, implicit, l=V_{DD}
D2 4 2; up
W 4 4_1; down=0.1, implicit, l=0V
W 1 2; right=0.5
W U1._VSS 4_2; down=0.6, implicit, l=0V
W U1._VDD 3_2; up=0.6, implicit, l=V_{DD}
#; help_lines=0.05
;draw_nodes=connections,label_nodes=false


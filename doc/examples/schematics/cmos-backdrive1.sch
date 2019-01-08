U1 inverter ._IN ._VSS ._OUT ._VDD; right
U2 inverter ._IN ._VSS ._OUT ._VDD; right
W U1._OUT U2._IN; right=1
W U1._VSS 4_2; down=0.3, implicit, l=0V
W U2._VSS 4_3; down=0.3, implicit, l=0V
W U1._VDD 3_2; up=0.3, implicit, l=V_{DD1}
W U2._VDD 3_3; up=0.3, implicit, l=V_{DD2}
#; help_lines=0.05
;draw_nodes=connections,label_nodes=false,thickness=2
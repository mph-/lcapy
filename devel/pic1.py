from schematic import Schematic

sch = Schematic()

sch.net_add('P1 1 0.1; v=V_1')
sch.net_add('R1 3 1; right')
sch.net_add('L1 2 3; right')
sch.net_add('C1 3 0; up')
sch.net_add('P2 2 0.2; v=V_2')

sch.draw()


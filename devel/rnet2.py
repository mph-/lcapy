from schematic import Schematic

sch = Schematic()

sch.net_add('R1 2 1')
sch.net_add('R2 3 2')
sch.net_add('R3 4 3')
sch.net_add('R4 6 1_1')
sch.net_add('R5 4_1 6')
sch.net_add('W1 1_1 1; up, size=0.5')
sch.net_add('W2 4_1 4; up, size=0.5')



sch.draw()


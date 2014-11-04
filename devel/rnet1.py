from schematic import Schematic

sch = Schematic()

sch.net_add('R1 2 1')
sch.net_add('R2 3 2')
sch.net_add('R3 4 3')
sch.net_add('R4 6 5')
sch.net_add('R5 7 6')
sch.net_add('W1 5 1; up')
sch.net_add('W2 7 4; up')



sch.draw()


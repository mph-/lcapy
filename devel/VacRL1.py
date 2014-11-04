from schematic import Schematic

sch = Schematic()

sch.net_add('Vac1 1 0.1; up')
sch.net_add('R1 2 1; right')
sch.net_add('L1 2 0; up')

sch.draw()


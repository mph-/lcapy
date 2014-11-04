from schematic import Schematic

sch = Schematic()

sch.add('P1 1 0.1; size=2')
sch.add('R1 3 1; right')
sch.add('L1 2 3; right')
sch.add('C1 3 4; up')
sch.add('L2 4 0; up')
sch.add('P2 2 0.2; size=2')

sch.draw()


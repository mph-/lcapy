from schematic import Schematic

sch = Schematic()

sch.add('Vac1 1 0.1; up')
sch.add('R1 2 1; right')
sch.add('L1 2 0; up')

sch.draw()


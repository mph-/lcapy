from schematic import Schematic

sch = Schematic()

sch.add('V1 1 0.1 ac; up')
sch.add('R1 2 1; right')
sch.add('L1 2 0; up')

sch.draw()


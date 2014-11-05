from schematic import Schematic

sch = Schematic()

sch.add('R1 2 1')
sch.add('R2 3 2')
sch.add('R3 4 3')
sch.add('R4 6 1_1')
sch.add('R5 4_1 6')
sch.add('W1 1_1 1; down, size=0.5')
sch.add('W2 4_1 4; down, size=0.5')



sch.draw()


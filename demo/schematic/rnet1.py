from lcapy import Schematic

sch = Schematic()

sch.add('R1 2 1')
sch.add('R2 3 2')
sch.add('R3 4 3')
sch.add('R4 6 5')
sch.add('R5 7 6')
sch.add('W1 5 1; down')
sch.add('W2 7 4; down')



sch.draw()


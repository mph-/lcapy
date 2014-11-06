from lcapy import Schematic

sch = Schematic()

sch.add('R1 1 2')
sch.add('R2 2 3')
sch.add('R3 3 4')
sch.add('R4 1_1 5')
sch.add('R5 5 4_1')
sch.add('W1 1_1 1; down, size=0.5')
sch.add('W2 4_1 4; down, size=0.5')



sch.draw(tex=True)


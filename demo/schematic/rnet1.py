from lcapy import Schematic

sch = Schematic()

sch.add('R1 1 2')
sch.add('R2 2 3')
sch.add('R3 3 4')
sch.add('R4 5 6')
sch.add('R5 6 7')
sch.add('W1 5 1; down')
sch.add('W2 7 4; down')



sch.draw(tex=True)


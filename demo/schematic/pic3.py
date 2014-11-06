from lcapy import Schematic

sch = Schematic()

sch.add('P1 1 0.1')
sch.add('R1 1 3; right')
sch.add('L1 3 2; right')
sch.add('C1 3 4; down')
sch.add('L2 4 0; down')
sch.add('P2 2 0.2')

sch.draw(tex=True)


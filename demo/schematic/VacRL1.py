from lcapy import Schematic

sch = Schematic()

sch.add('V1 1 0.1 ac; down')
sch.add('R1 1 2; right')
sch.add('L1 2 0; down')

sch.draw(tex=True)


from lcapy import Schematic

sch = Schematic()

sch.add('V1 1 0.1; down')
sch.add('R1 1 2; right, i=I_1, v=V_{R_1}, size=1.5')
sch.add('L1 2 0; down, l=L_{23}, i=I_1, v=V_{L_1}, size=1.5')
sch.add('W  0.1 0; right')

sch.draw(tex=True)


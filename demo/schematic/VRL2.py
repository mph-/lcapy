from lcapy import Schematic

sch = Schematic()

sch.add('V1 1 0; down')
sch.add('R1 1 2; left, i=I_1, v=V_{R_1}')
sch.add('R2 1 3; right, i=I_2, v=V_{R_2}')
sch.add('L1 2 0.1; down, i=I_1, v=V_{L_1}')
sch.add('L2 3 0.3; down, i=I_1, v=V_{L_2}')
sch.add('W 0 0.3; right')
sch.add('W 0 0.1; left')

sch.draw(scale=3)


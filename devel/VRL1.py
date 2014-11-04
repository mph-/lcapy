from schematic import Schematic

sch = Schematic()

sch.add('V1 1 0.1; up')
sch.add('R1 2 1; right, i_=I_1, v=V_x')
sch.add('L1 2 0; up, l=L_{23}')

sch.draw()


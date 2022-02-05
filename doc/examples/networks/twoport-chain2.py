from lcapy import Shunt, Series, R, V

tp = Shunt(R('R1')).chain(Series(R('R2') + V('V1')))

tp.draw(__file__.replace('.py', '.png'))

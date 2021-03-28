from lcapy import Circuit, f, oo
from numpy import logspace

a = Circuit('opamp-piezo-amplifier1.sch')

b = a.subs({'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':100e-9, 'Vn':2e-9, 'Inp':5e-15, 'Inn':5e-15})

Vo = b.Po.V.n(f).limit('A', oo).simplify()

vf = logspace(0, 5, 201)

Vo.plot(vf, log_frequency=True, yscale=1e9, ylabel='Voltage noise ASD (nV/rtHz)')

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')




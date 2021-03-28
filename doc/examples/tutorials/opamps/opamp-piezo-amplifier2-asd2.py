from lcapy import Circuit, f, oo
from numpy import logspace

a = Circuit('opamp-piezo-amplifier1.sch')

Vo = a.Po.V.n(f).limit('A', oo)

defs1 = {'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':100e-6, 'Vn':2e-9, 'Inp':5e-15, 'Inn':5e-15}
defs2 = {'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':10e-6, 'Vn':2e-9, 'Inp':5e-15, 'Inn':5e-15}
defs3 = {'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':1e-6, 'Vn':2e-9, 'Inp':5e-15, 'Inn':5e-15}

Vo1 = Vo.subs(defs1)
Vo2 = Vo.subs(defs2)
Vo3 = Vo.subs(defs3)

vf = logspace(0, 5, 201)

ax = Vo1.plot(vf, log_frequency=True, yscale=1e9, label='C=100 uF' )
ax = Vo2.plot(vf, log_frequency=True, yscale=1e9, label='C=10 uF', axes=ax)
Vo3.plot(vf, log_frequency=True, yscale=1e9, label='C=1 uF', ylabel='Voltage noise ASD (nV/rtHz)', axes=ax)
ax.legend()

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')




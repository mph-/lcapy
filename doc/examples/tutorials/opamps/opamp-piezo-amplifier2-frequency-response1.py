from lcapy import Circuit, f, oo
from numpy import logspace

a = Circuit('opamp-piezo-amplifier1.sch')

H = a.transfer('Cs', 'Po')(f).limit('A', oo)

defs1 = {'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':100e-6, 'Vn':2e-9, 'Inp':5e-15, 'Inn':5e-15}
defs2 = {'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':10e-6, 'Vn':2e-9, 'Inp':5e-15, 'Inn':5e-15}
defs3 = {'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':1e-6, 'Vn':2e-9, 'Inp':5e-15, 'Inn':5e-15}

H1 = H.subs(defs1)
H2 = H.subs(defs2)
H3 = H.subs(defs3)

vf = logspace(0, 5, 201)

ax = abs(H1).plot(vf, log_frequency=True, label='C=100 uF' )
ax = abs(H2).plot(vf, log_frequency=True, label='C=10 uF', axes=ax)
abs(H3).plot(vf, log_frequency=True, label='C=1 uF', ylabel='Gain (dB)', axes=ax)
ax.legend()

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')




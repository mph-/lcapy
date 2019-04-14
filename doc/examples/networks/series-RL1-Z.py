from lcapy import *
from numpy import logspace
from matplotlib.pyplot import figure, savefig

N = R(10) + L(1e-2)

vf = logspace(0, 5, 400)
Z = N.Z(f).evaluate(vf)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(vf, abs(Z), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)

savefig('series-RL1-Z.png')


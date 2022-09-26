from lcapy import *
from numpy import logspace
from matplotlib.pyplot import figure, savefig

N = R(10) | C(1e-4) | L(1e-3)

vf = logspace(0, 5, 400)
Z = N.Z(j2pif).evaluate(vf)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(vf, abs(Z), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)


savefig('parallel-RLC3-Z.png')

from lcapy import *
from numpy import logspace
from matplotlib.pyplot import figure, savefig, show

N = R(10) + C(1e-4)

vf = logspace(0, 5, 400)
Z = N.Z.frequency_response(vf)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(vf, abs(Z), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)
show()

savefig('series-RC1-Z.png')

from mcircuit import *
from numpy import logspace
from matplotlib.pyplot import figure, savefig, show

N = V(10) + R(10) + C(1e-4)

f = logspace(0, 5, 400)
Isc = N.Isc.frequency_response(f)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(Isc), linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A/Hz)')
ax.grid(True)
show()

savefig('series-VRC1-Isc-frequency.png')

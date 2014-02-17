from mcircuit import *
from numpy import linspace

N = V(20) + R(5) + C(10)

t = linspace(0, 100, 400)
isc = N.Isc.transient_response(t)

from matplotlib.pyplot import figure, savefig, show
fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, isc, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)
show()

savefig('series-RLC1-isc.png')

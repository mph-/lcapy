from mcircuit import *
from numpy import logspace
from matplotlib.pyplot import figure, savefig, show

N = V(20) + R(5) + C(10)

t = np.linspace(0, 100, 1000)
isc = N.Isc.transient_response(t)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, isc, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)
show()

savefig('series-VRC1-isc.png')

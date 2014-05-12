from lcapy import *
from numpy import logspace
from matplotlib.pyplot import figure, savefig, show

N = Vdc(20) + R(10) + C(1e-4)

t = np.linspace(0, 0.01, 1000)
isc = N.Isc.transient_response(t)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, isc, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)
show()

savefig('series-VRC1-isc.png')

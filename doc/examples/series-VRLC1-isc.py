from mcircuit import V, R, L, C
from matplotlib.pyplot import figure, savefig, show
import numpy as np

a = V(10) + R(0.1) + C(0.4) + L(0.2)

a.Isc.pprint()

t = np.linspace(0, 10, 1000)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, a.Isc.transient_response(t), linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)

savefig('series-VRLC1-isc.png')

show()



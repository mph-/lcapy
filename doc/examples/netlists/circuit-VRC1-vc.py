from lcapy import Circuit

cct = Circuit()
cct.add('V 1 0 step 20')
cct.add('R 1 2 10')
cct.add('C 2 0 1e-4')

# Determine transient response at node 2 evaluated at specified times.
import numpy as np
t = np.linspace(0, 0.01, 1000)
vc = cct[2].v.evaluate(t)

from matplotlib.pyplot import figure, savefig
fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, vc, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Capacitor voltage (V)')
ax.grid(True)

savefig('circuit-VRC1-vc.png')

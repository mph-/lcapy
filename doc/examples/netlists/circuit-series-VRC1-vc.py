from lcapy import *
import numpy as np
from matplotlib.pyplot import subplots, savefig

t = np.linspace(0, 0.01, 1000)

cct = Circuit()

cct.add('V1 1 0 step 20')
cct.add('R1 1 2 10')
cct.add('C1 2 0 1e-4')

vc = cct.C1.v.evaluate(t)

fig, ax = subplots(1)
ax.plot(t, vc, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Capacitor voltage (V)')
ax.grid(True)

savefig('circuit-series-VRC1-vc.png')

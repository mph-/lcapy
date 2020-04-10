from lcapy import Circuit

cct = Circuit("""
V 1 0 step 20
R 1 2 10
C 2 0 1e-4
""")

import numpy as np
t = np.linspace(0, 0.01, 1000)
vc = cct.C.v.evaluate(t)

from matplotlib.pyplot import subplots, savefig
fig, ax = subplots(1)
ax.plot(t, vc, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Capacitor voltage (V)')
ax.grid(True)

savefig('circuit-VRC1-vc.png')

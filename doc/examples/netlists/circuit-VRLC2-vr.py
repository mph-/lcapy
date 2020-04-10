from lcapy import Circuit
cct = Circuit("""
V 1 0 step 10; down
L 1 2 1e-3; right, size=1.2
C 2 3 1e-4; right, size=1.2
R 3 0_1 1; down
W 0 0_1; right
""")

import numpy as np
t = np.linspace(0, 0.01, 1000)
vr = cct.R.v.evaluate(t)

from matplotlib.pyplot import subplots, savefig
fig, ax = subplots(1)
ax.plot(t, vr, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Resistor voltage (V)')
ax.grid(True)

savefig('circuit-VRLC2-vr.png')

from lcapy import Circuit
cct = Circuit()

cct.add('V 1 0 step 10')
cct.add('L 1 2 1e-3')
cct.add('C 2 3 1e-4')
cct.add('R 3 0 10')

import numpy as np
t = np.linspace(0, 0.01, 1000)
vr = cct.R.v(t)

from matplotlib.pyplot import figure, savefig, show
fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, vr, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Resistor voltage (V)')
ax.grid(True)
show()

savefig('circuit-VRLC1-vr.png')

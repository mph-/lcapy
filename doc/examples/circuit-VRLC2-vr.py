from lcapy import Circuit
import numpy as np
from matplotlib.pyplot import figure, savefig, show

t = np.linspace(0, 0.01, 1000)

cct = Circuit()

cct.add('V1 1 0 dc 10')
cct.add('L1 1 2 1e-3')
cct.add('C1 2 3 1e-4')
cct.add('R1 3 0 1')

Vr = cct[3].V
vr = Vr.transient_response(t)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, vr, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Resistor voltage (V)')
ax.grid(True)
show()

savefig('circuit-VRLC2-vr.png')

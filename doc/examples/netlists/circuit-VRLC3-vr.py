from lcapy import Circuit
import numpy as np
from matplotlib.pyplot import subplots, savefig

t = np.linspace(0, 0.01, 1000)

cct = Circuit()

cct.add('V1 1 0 step 10')
cct.add('L1 1 2 0.1')
cct.add('C1 2 3 1e-4')
cct.add('R1 3 0 100')

Vr = cct[3].V
vr = Vr.transient_response(t)

fig, ax = subplots(1)
ax.plot(t, vr, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Resistor voltage (V)')
ax.grid(True)


savefig('circuit-VRLC3-vr.png')

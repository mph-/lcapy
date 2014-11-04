from lcapy import *
import numpy as np
from matplotlib.pyplot import figure, savefig, show

t = np.linspace(0, 0.01, 1000)

cct = Circuit('Series VRC1')

cct.add('V1 1 0 dc 20')
cct.add('R1 1 2 10')
cct.add('C1 2 0 1e-4')


Vc = cct.V[2]
vc = Vc.transient_response(t)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, vc, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Capacitor voltage (V)')
ax.grid(True)
show()

savefig('circuit-series-VRC1-vc.png')

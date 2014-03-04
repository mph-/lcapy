from mcircuit import V, R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

Vbody = 5e3
Cbody = 100e-12
Rbody = 1.5e3
Cdev = 5e-12


a1 = V(Vbody) + C(Cbody) + R(Rbody)
b1 = a1.load(C(Cdev))

t = np.linspace(0, 50e-9, 1000)

fig = figure()
ax = fig.add_subplot(111)
# Voltage across R
ax.plot(t * 1e9, b1.V.impulse_response(t) / 1e3, linewidth=2)
ax.legend()
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Voltage (kV)')
ax.grid(True)

vdev = b1.V.impulse_response(t)
idev = b1.I.impulse_response(t)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(t * 1e9, vdev * idev / 1e3, linewidth=2)
ax.legend()
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Power (kW)')
ax.grid(True)


show()

from lcapy import Vdc, R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

# Human body model from MIL-STD-883, Method 3015.8, 
Cbody = 100e-12
Rbody = 1.5e3

# Voltage on body.
Vbody = 5e3

# Device input capacitance.
Cdev = 5e-12


a1 = Vdc(Vbody) + C(Cbody) + R(Rbody)
b1 = a1.load(C(Cdev))

t = np.linspace(0, 50e-9, 1000)

fig = figure()
ax = fig.add_subplot(111)
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

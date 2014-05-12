from lcapy import Vdc, R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

a = (Vdc(5) + L(10)) | R(5)
b = Vdc(5) + L(10) + R(5)

t = np.linspace(0, 10, 1000)

fig = figure()
ax = fig.add_subplot(111)
# Voltage across R
ax.plot(t, a.V.impulse_response(t), linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Voltage (V)')
ax.grid(True)


fig = figure()
ax = fig.add_subplot(111)
# Current through R
ax.plot(t, b.I.impulse_response(t), linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)

show()

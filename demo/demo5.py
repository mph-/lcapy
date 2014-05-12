from lcapy import Vdc, R, L, C, LSection, Shunt
import numpy as np
from matplotlib.pyplot import figure, savefig, show

a1 = LSection(L(10), C(1, 5))
b1 = a1.prepend(Shunt(Vdc(5))).load(R(5))

a2 = LSection(L(10), C(1e-10, 5))
b2 = a2.prepend(Shunt(Vdc(5))).load(R(5))

a1 = (Vdc(5) + L(10))
a2 = (Vdc(5) + L(10)) | C(1, 5)
b1 = a1.load(R(5))
b2 = a2.load(R(5))

t = np.linspace(0, 10, 1000)

fig = figure()
ax = fig.add_subplot(111)
# Voltage across R
ax.plot(t, b1.V.impulse_response(t), linewidth=2, label='without C')
ax.plot(t, b2.V.impulse_response(t), linewidth=2, label='with C')
ax.legend()
ax.set_xlabel('Time (s)')
ax.set_ylabel('Voltage (V)')
ax.grid(True)


fig = figure()
ax = fig.add_subplot(111)
# Current through R
ax.plot(t, b1.I.impulse_response(t), linewidth=2, label='without C')
ax.plot(t, b2.I.impulse_response(t), linewidth=2, label='with C')
ax.legend()
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)

show()


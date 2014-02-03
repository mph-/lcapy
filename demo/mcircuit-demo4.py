from msignal.mcircuit import V, R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

a1 = (V(5) + L(10))
a2 = (V(5) + L(10)) | C(1, 5)
b1 = a1.load(R(5))
b2 = a2.load(R(5))

t = np.linspace(0, 10, 1000)

fig = figure()
ax = fig.add_subplot(111)
# Voltage across R
ax.plot(t, b1.V.impulseresponse(t), linewidth=2, label='without C')
ax.plot(t, b2.V.impulseresponse(t), linewidth=2, label='with C')
ax.legend()
ax.set_xlabel('Time (s)')
ax.set_ylabel('Voltage (V)')
ax.grid(True)


fig = figure()
ax = fig.add_subplot(111)
# Current through R
ax.plot(t, b1.I.impulseresponse(t), linewidth=2, label='without C')
ax.plot(t, b2.I.impulseresponse(t), linewidth=2, label='with C')
ax.legend()
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)

show()

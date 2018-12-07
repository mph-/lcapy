from lcapy import Vstep, R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

# This is posed as an initial value problem so cannot
# determine result for t < 0.
a1 = Vstep(5) + L(10, 0)
a2 = a1 | C(1, 5)
b1 = a1 | R(5)
b2 = a2 | R(5)

t = np.linspace(-1, 10, 1000)

fig = figure()
ax = fig.add_subplot(111)
# Open-circuit voltage across R
ax.plot(t, b1.v.evaluate(t), linewidth=2, label='without C')
ax.plot(t, b2.v.evaluate(t), linewidth=2, label='with C')
ax.legend()
ax.set_xlabel('Time (s)')
ax.set_ylabel('Voltage (V)')
ax.grid(True)


fig = figure()
ax = fig.add_subplot(111)
# Short-circuit current through R
ax.plot(t, b1.isc.evaluate(t), linewidth=2, label='without C')
ax.plot(t, b2.isc.evaluate(t), linewidth=2, label='with C')
ax.legend()
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current (A)')
ax.grid(True)

show()

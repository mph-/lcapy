from lcapy import Vstep, R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

a = (Vstep(5) + L(10)) | R(5)
b = Vstep(5) + L(10) + R(5)

t = np.linspace(-1, 10, 1000)

fig = figure()
ax = fig.add_subplot(111)
# Open-circuit voltage across R
ax.plot(t, a.v.evaluate(t), linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Voltage (V)')
ax.grid(True)


fig = figure()
ax = fig.add_subplot(111)
# Short-circuit current through R
ax.plot(t, b.isc.evaluate(t), linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Short-circuit current (A)')
ax.grid(True)

show()

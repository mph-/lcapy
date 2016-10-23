from lcapy import Istep, R, L, C
from matplotlib.pyplot import figure, savefig, show
import numpy as np

a = Istep(10) | R(0.1) | C(0.4) | L(0.2)

a.Voc.pprint()

t = np.linspace(0, 10, 1000)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(t, a.Voc.transient_response(t), linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Voltage (V)')
ax.grid(True)

savefig('parallel-IRLC1-voc.png')

show()



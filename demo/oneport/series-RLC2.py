from lcapy import Vstep, R, L, C
from matplotlib.pyplot import figure, savefig, show
import numpy as np

a = R(0.1) + C(0.4) + L(0.2)

a.Z.pprint()

f = np.linspace(0, 1000, 1000)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(f, abs(a.Z.frequency_response(f)), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)


show()

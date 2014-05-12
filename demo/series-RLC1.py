from lcapy import Vdc, R, L, C
from matplotlib.pyplot import figure, savefig, show
import numpy as np

a = R(10) + C(1e-4) + L(1e-3)

a.Z.canonical().pprint()
a.Z.ZPK().pprint()

f = np.logspace(0, 5, 1000)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(a.Z.frequency_response(f)), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)


show()


b = R('R') + C('C') + L('L')

b.Z.canonical().pprint()
b.Z.ZPK().pprint()

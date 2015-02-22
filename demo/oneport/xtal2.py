from lcapy import R, L, C, Xtal, pprint
import numpy as np
from matplotlib.pyplot import figure, savefig, show

C_0 = 4e-12
f_1 = 10e6
C_1 = 8e-15
L_1 = 1 / ((2 * np.pi * f_1)**2 * C_1)
R_1 = 20

xtal = Xtal(C_0, R_1, L_1, C_1)

pprint(xtal)
pprint(xtal.expand())
pprint(xtal.Z)

f = np.logspace(6, 8, 1000)

fig = figure()
ax = fig.add_subplot(111)
Zf = xtal.Z.frequency_response(f)
ax.loglog(f, abs(Zf))
ax.grid(True)

show()

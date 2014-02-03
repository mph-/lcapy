from __future__ import division
from mcircuit import R, L, C, Xtal
import numpy as np
from matplotlib.pyplot import figure, savefig, show

C_0 = 4e-12
f_1 = 10e6
C_1 = 8e-15
L_1 = 1 / ((2 * np.pi * f_1)**2 * C_1)
R_1 = 20

xtal = Xtal(C_0, R_1, L_1, C_1)

ZCb = C(20e-12)

H = ZCb.divider(xtal).H

f = np.logspace(6, 8, 2000)
Hf = H.freqresponse(f)
Zfxtal = xtal.Z.freqresponse(f)

fig = figure()
ax = fig.add_subplot(111)
ax.semilogx(f, np.angle(Hf) / np.pi * 180)
ax.grid(True)

fig = figure()
ax = fig.add_subplot(111)
ax.semilogx(f, abs(Hf))
ax.grid(True)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(Zfxtal))
ax.grid(True)


show()

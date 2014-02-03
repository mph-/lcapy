from __future__ import division
from msignal.mcircuit import R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

ZC_0 = C(10e-12)

f_1 = 10e6
C_1 = 8e-15
L_1 = 1 / ((2 * np.pi * f_1)**2 * C_1)

ZL_1 = L(L_1)
ZR_1 = R(10 * 2)
ZC_1 = C(C_1)

Zxtal = (ZL_1 + ZR_1 + ZC_1).parallel(ZC_0)

ZC_b = C(20e-12)

ZR_f = R(1e6)

H = ZC_b.divider(Zxtal.parallel(ZR_f))

f = np.logspace(6, 8, 2000)

fig = figure()
ax = fig.add_subplot(111)
Hf = H.freqresponse(f)
ax.semilogx(f, np.angle(Hf) / np.pi * 180)
ax.grid(True)


fig = figure()
ax = fig.add_subplot(111)
Hf = H.freqresponse(f)
ax.semilogx(f, abs(Hf))
ax.grid(True)


show()

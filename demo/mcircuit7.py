from __future__ import division
from mcircuit import R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

ZC = C(40e-12)
ZL = L(0.06e-4)
ZR = R(10)

H = ZC.divider(ZR + ZL)

f = np.logspace(6, 8, 2000)

fig = figure()
ax = fig.add_subplot(111)
Hf = H.frequency_response(f)
ax.semilogx(f, np.angle(Hf) / np.pi * 180)
ax.grid(True)


fig = figure()
ax = fig.add_subplot(111)
Hf = H.frequency_response(f)
ax.semilogx(f, abs(Hf))
ax.grid(True)


show()

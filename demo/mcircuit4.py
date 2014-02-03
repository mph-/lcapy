#import msignal.mcircuit as mc
from __future__ import division
from msignal.mcircuit import R, L, C
import numpy as np
from matplotlib.pyplot import figure, savefig, show

ZC_0 = C(4e-12)

f_1 = 20e6
C_1 = 8e-15
L_1 = 1 / ((2 * np.pi * f_1)**2 * C_1)

ZL_1 = L(L_1)
ZR_1 = R(10 * 2)
ZC_1 = C(C_1)

Zxtal = (ZL_1 + ZR_1 + ZC_1).parallel(ZC_0)

ZC_a = C(20e-9)
ZC_b = C(20e-9)
ZC_i = C(5e-12)
ZC_o = C(5e-12)
ZR_o = R(20)
ZR_lim = R(1000)
ZR_f = R(1e6)

f = np.logspace(6, 8, 2000)

Ztot = ZR_f + ZR_lim + Zxtal

# Use delta to wye transform.
H = ZC_i.parallel(ZC_a).chain_divider(ZR_f * Zxtal / Ztot, ZR_lim * Zxtal / Ztot + ZC_b, ZR_f * ZR_lim / Ztot)

fig = figure()
ax = fig.add_subplot(111)
Hf = H.freqresponse(f)
ax.semilogx(f, np.angle(Hf) / np.pi * 180)
ax.grid(True)


show()

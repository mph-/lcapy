# This is the same as loudspeaker1.py but using the topology changing
# mechanical model.

from lcapy import V, R, L, C, IdealTransformer, Shunt, Series, TSection, IdealGyrator
import numpy as np
from matplotlib.pyplot import figure, savefig, show

Re = 5.4
Le = 3.3e-3
fs = 21.5
Qms = 2.37
Mms = 54.8e-3
Cms = 1.0e-3
S = 310e-4
Bl = 10.0

omegas = 2 * np.pi * fs
Rms = Qms / (omegas * Cms)

Rs = Bl**2 / Rms
Ls = Bl**2 * Cms
Cs = Mms / (Bl ** 2)

print Rms
print Rs
print Ls
print Cs

a = Series(R(Re) + L(Le)).chain(IdealTransformer(1 / Bl))
b = a.chain(Shunt(R(1 / Rms) | C(Mms) | L(Cms)))

f = np.logspace(0, 5, 2000)
Zin = b.open_circuit(2).Z.frequency_response(f)

fig = figure()
ax = fig.add_subplot(111)
ax.semilogx(f, Zin, linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)

show()





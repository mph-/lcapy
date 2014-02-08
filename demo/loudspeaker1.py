# This is the same as loudspeaker1.py but using the topology consistent
# mechanical model.

from mcircuit import V, R, L, C, IdealTransformer, Shunt, Series, TSection, IdealGyrator
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

a = Series(R(Re) + L(Le)).chain(IdealGyrator(Bl))
b = a.chain(Series(R(Rms) + L(Mms) + C(Cms)))

f = np.logspace(0, np.log10(20e3), 2000)
Zin = b.shortcircuit(2).Z.frequency_response(f)

fig = figure()
ax = fig.add_subplot(111)
#ax.semilogx(f, abs(Zin), linewidth=2)
ax.semilogx(f, Zin.real, linewidth=2)
#ax.semilogx(f, Zin.imag, linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Resistance (ohms)')
ax.grid(True)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(Zin), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance (ohms)')
ax.grid(True)

show()





from mcircuit import R, L, C, LSection, Shunt
import numpy as np
from matplotlib.pyplot import figure, savefig, show

# Gibson PAF pickup
L0 = 4.4

# 42 AWG 1.6 ohms/foot
# 43 AWG 2.1 ohms/foot

# 5000 turns, 42 AWG, 1.6 ohms/foot

f0 = 7.715e3
R0 = 319e3
# The impedance at resonance will be slightly higher than R0.

omega0 = 2 * np.pi * f0

C0 = 1 / (omega0**2 * L0)

a = LSection(R(R0) + L(L0), C(C0))

R1 = 500e3
R1 = 1e3
C1 = 22e-9

a = a.chain(Shunt(R(R1) + C(C1)))

f = np.logspace(1, 5, 1000)

H = a.Vgain12
Hf = H.frequency_response(f)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(Hf), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.grid(True)

b = a.short_circuit(1)

Zoc = b.Zoc
Zocf = Zoc.frequency_response(f)


fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(Zocf), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.grid(True)

c = (R(R0) + L(L0)) | C(C0)


Zoc = c.Zoc
Zocf = Zoc.frequency_response(f)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(Zocf), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.grid(True)


show()




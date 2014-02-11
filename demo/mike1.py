from msignal.mcircuit import V, R, L, C, IdealTransformer, Shunt, Series, TSection
import numpy as np
from matplotlib.pyplot import figure, savefig, show

# The following parameters are from the SD02 datasheet for a 60 mm
# spike.  It appears the sensor relates to the aluminium part.  The
# total includes the steel spike.

# Charge sensivity (sensor) typ at 160 Hz
QS = 1.35e-12

# Force sensivity (total)   typ at 160 Hz
FS = 23e-12

# Capacitance (without cable)  +/- 5%
C0 = 0.26e-9

# Cable capacitance (for 1 m)
Cc = 100e-12 * 1

# Mass (sensor) (kg)
ms = 22e-3

# Mass (total) (kg)
mt = 62e-3

# Resonant frequency (Hz)  (this is stated for the sensor and NA for total)
f0 = 23e3

# The following parameters are measured or guessed.

# Electromechanical factor (N/V)  guessed for PZT
alpha = 0.02
# The following gives a better match
alpha = 0.027

# Mechanical resistance (kg/s) from Simon Woods?
# This will only affect the damping zeta
Rm = alpha**2 / 16.4e-6

# Mass of brass load from Simon Woods
m1 = 1.8e-3

# Mass of crystal from Simon Woods
mx = 3.0e-3

# Mass of body (kg)
m2 = mt - m1 - mx

# Effective masses (kg)
m1p = m1 + mx / 2
m2p = m2 + mx / 2

# Piezo compliance (m/N)
Sm = 1 / (m1p * (2 * np.pi * f0)**2)
Sm = 1 / ((2 * np.pi * f0)**2 * m1p - alpha**2 / C0)

# Gravitational acceleration (m/s^2)
g = 9.81

# Sensitivity V/m/s^2
H = m1p * Sm * alpha / C0

print('Piezo conversion factor %.3f N/V' % alpha)
print('Sensitivity %.1f mV/m/s^2 (%.3f mV/g)' % (H * 1e3, H / g * 1e3))

# Sensitivity V/N
G = H / m2p

print('Sensitivity %.1f mV/N' % (G * 1e3))

print('Force sensitivity %.1f pC/N cf. %.1f' % (G * C0 * 1e12, FS * 1e12))

print('Charge sensitivity %.1f pC/m/s^2 cf. %.1f' % (H * C0 * 1e12, QS * 1e12))


print("Masses m1' %.1f g  m2' %.1f g" % (m1p * 1e3, m2p * 1e3))
print('Compliance Sm %.1f nm/N' % (Sm * 1e9))
print('Resistance Rm %.1f kg/s' % Rm)

omega0 = np.sqrt ((1 / Sm + alpha**2 / C0) / m1p)
f0 = omega0 / (2 * np.pi)
zeta = Rm / m1p / (2 * omega0)
K = alpha  / C0

omega1 = omega0 * np.sqrt (1 - zeta**2)
f1 = omega1 / (2 * np.pi)

print('Resonant frequency f0 %.2f kHz' % (f0 * 1e-3))
print('Damped resonant frequency f1 %.2f kHz' % (f1 * 1e-3))
print('Damping factor %.3f' % zeta)

me = m1p * m2p / (m1p + m2p)
L1 = me / alpha**2
C1 = Sm * alpha**2
R1 = Rm / alpha**2

print('C0 %.1f pF' % (C0 * 1e12))
print('C1 %.1f pF' % (C1 * 1e12))
print('L1 %.1f H' % (L1 * 1))
print('R1 %.1f kohm' % (R1 * 1e-3))


# Let's say we measure 2 V with a gain of 80 and with the cable connected.
# This corresponds to a sensor voltage of 2 / 80 = 25 mV corresponding
# to an input force of 0.4 N (6.7 m/s^2).
Vo = 2.0 / 80
Fi = Vo / G
Ai = Fi / m2p
Qo = (C0 + Cc) * Vo

print('Developed charge %.1f pC for voltage %.1f mV' % (Qo * 1e12, Vo * 1e3))

# Wood compliance approximating radiation impedance 
Smw = 725e-9

a = TSection(L(m2p) | C(Smw), L(m1p), C(Sm) + R(Rm)).chain(IdealTransformer(1 / alpha)).chain(Shunt(C(C0)))

f = np.logspace(1, 5, 2000)
Av = a.Vgain12.freqresponse(f)
AvdB = 20 * np.log10(abs(Av))

fig = figure()
ax = fig.add_subplot(111)
ax.semilogx(f, AvdB, linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Response (dB)')
ax.grid(True)

show()




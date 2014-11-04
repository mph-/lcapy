from lcapy import pprint, Circuit
import numpy as np
from matplotlib.pyplot import figure, savefig, show
from msignal import Msignal

cct = Circuit('Transducer')

fr = 36e3
C1 = 1e-6
L1 = 1 / ((2 * np.pi * fr)**2 * C1)
R1 = 0.1

cct.add('Vi 1 0 dc') 
cct.add('C1 2 1', C1) 
cct.add('L1 3 2', L1) 
cct.add('R1 4 3', R1) 
cct.add('Rr 4 0 0.01') 



H = cct.V[4] / cct.V[1]
pprint(H)

f = np.logspace(np.log10(1e3), np.log10(1e6), 2000)
Hf = H.frequency_response(f)

fig = figure()
ax = fig.add_subplot(111)
ax.loglog(f, abs(Hf), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Gain')
ax.grid(True)


# f = np.linspace(1e3, 100e3, 2000)
# Hf = H.frequency_response(f)
#
# fig = figure()
# ax = fig.add_subplot(111)
# ax.semilogy(f, abs(Hf), linewidth=2)
# ax.set_xlabel('Frequency (Hz)')
# ax.set_ylabel('Gain')
# ax.grid(True)


Tp = 10e-3
f0 = 37.5e3

t = np.linspace(0, 1, 100000)
p = Msignal(np.sin(2 * np.pi * f0 * t) * (t < Tp), t)

P = p.fft()


fig = figure()
ax = fig.add_subplot(111)
ax.loglog(P.f, abs(P), linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Foo')
ax.grid(True)

pprint(H.transient_response())

h = H.transient_response(t)

e = Msignal(H.response(p, t)[0 : len(t)], t)

fig = figure()
ax = fig.add_subplot(111)
ax.plot(e.t, e, linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Amplitude')
ax.grid(True)

E = e.fft()

AHf = abs(Hf) / max(abs(Hf))

fig = figure()
ax = fig.add_subplot(111)
ax.plot(P.f, abs(P.norm()), linewidth=2)
ax.plot(E.f, abs(E.norm()), linewidth=2)
ax.plot(f, AHf, linewidth=2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Foo')
ax.set_xlim((32e3, 40e3))

ax.grid(True)


show()

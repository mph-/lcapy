from lcapy import s, j, pi, f, transfer, j2pif
from matplotlib.pyplot import savefig
from numpy import logspace

# Note, this has a marginally stable impulse response since it has a
# pole at s = 0.
H = transfer((s - 2) * (s + 3) / (s * (s - 2 * j) * (s + 2 * j)))

fv = logspace(-1, 3, 400)
H(j2pif).dB.plot(fv, log_scale=True)

savefig('tf1-bode-plot.png')

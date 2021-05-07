from lcapy import s, j, pi, f, transfer
from matplotlib.pyplot import savefig
from numpy import logspace

H = transfer((s - 2) * (s + 3) / (s * (s - 2 * j) * (s + 2 * j)))

fv = logspace(-1, 3, 400)
H(f).dB.plot(fv, log_scale=True)

savefig('tf1-bode-plot.png')

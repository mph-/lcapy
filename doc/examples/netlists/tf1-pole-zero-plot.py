from lcapy import s, j, transfer
from matplotlib.pyplot import savefig

H = transfer((s - 2) * (s + 3) / (s * (s - 2 * j) * (s + 2 * j)))
H.plot()

savefig('tf1-pole-zero-plot.png')

from lcapy import s, j, Hs
from matplotlib.pyplot import savefig, show

H = Hs((s - 2) * (s + 3) / (s * (s - 2 * j) * (s + 2 * j)))
H.plot()

show()
savefig('tf1-pole-zero-plot.png')

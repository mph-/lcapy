from lcapy import delta
from lcapy.discretetime import n, z
from matplotlib.pyplot import savefig

x = delta(n) + delta(n - 2)
X = x(z)
X.plot()

savefig('dt1-pole-zero-plot1.png')

from lcapy import ui
from lcapy.discretetime import n, z
from matplotlib.pyplot import savefig

x = ui(n) + ui(n - 2)
X = x(z)
X.plot()

savefig('dt1-pole-zero-plot1.png')

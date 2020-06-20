from lcapy import delta
from lcapy.discretetime import n
from matplotlib.pyplot import savefig

x = delta(n) + delta(n - 2)
x.plot(figsize=(6, 2))

savefig('dt1-plot1.png')

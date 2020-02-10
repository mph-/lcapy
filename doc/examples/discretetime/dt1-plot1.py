from lcapy import ui
from lcapy.discretetime import n
from matplotlib.pyplot import savefig

x = ui(n) + ui(n - 2)
x.plot(figsize=(6, 2))

savefig('dt1-plot1.png')

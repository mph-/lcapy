from lcapy import ui
from lcapy.discretetime import n, dt
from matplotlib.pyplot import savefig

x = ui(n) + ui(n - 2)
x.DTFT().subs(dt, 1).plot(figsize=(6, 3))

savefig('dt1-DTFT-plot1.png', bbox_inches='tight')

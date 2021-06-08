from lcapy import delta
from lcapy.discretetime import n, dt
from matplotlib.pyplot import savefig

x = delta(n) + delta(n - 2)
abs(x.DTFT().subs(dt, 1)).plot(norm=True, figsize=(6, 3))

savefig('dt1-DTFT-plot1.png', bbox_inches='tight')

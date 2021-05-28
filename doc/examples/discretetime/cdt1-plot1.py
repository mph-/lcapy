from lcapy import j, n, exp
from matplotlib.pyplot import savefig

x = 0.9**n * exp(j * n * 0.5)
x.plot((1, 10), figsize=(6, 6), polar=True)

savefig('cdt1-plot1.png')

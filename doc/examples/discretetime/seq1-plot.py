from lcapy import seq
from matplotlib.pyplot import savefig

x = seq((1, 2, 1, 0))
x.plot(figsize=(6, 2), unfilled=True)

savefig('seq1-plot.png')

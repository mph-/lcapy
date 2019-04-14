from lcapy import Circuit

cct = Circuit("""
V 1 0 step 20
R 1 2 10
C 2 0 1e-4
""")

vc = cct.C.v

from matplotlib.pyplot import figure, savefig
import numpy as np

t = np.linspace(0, 0.01, 1000)

fig = figure()
ax = fig.add_subplot(111)
vc.plot(t, axes=ax)

savefig('circuit-VRC2-vc.png')

from lcapy import Circuit

cct = Circuit()
cct.add('V 1 0 step 20')
cct.add('R 1 2 10')
cct.add('C 2 0 1e-4')

# Determine transient response at node 2.
vc = cct[2].v

from matplotlib.pyplot import figure, savefig
import numpy as np

t = np.linspace(0, 0.01, 1000)

fig = figure()
ax = fig.add_subplot(111)
vc.plot(t, axes=ax)

savefig('circuit-VRC2-vc.png')

from msignal.mcircuit import Opamp, R, C, Series
import numpy as np
from matplotlib.pyplot import figure, savefig, show

# Create simple integrator
Cf = 1e-6
Ri = 50

a = Opamp()

# Connect V+ to ground.
b = a.shortcircuit(1)

# Add feedback capacitor.
c = b.bridge(C(Cf))

# Add input resistor to V-.
d = c.prepend(Series(R(Ri)))

print d.Vtransfer

f = np.logspace(1, 8, 1000)

fig = figure()
ax = fig.add_subplot(111)
Zf = d.Vtransfer.freqresponse(f)
ax.loglog(f, abs(Zf))
ax.grid(True)

show()

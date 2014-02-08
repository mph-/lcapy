from mcircuit import Opamp, R, C, Series
import numpy as np
from matplotlib.pyplot import figure, savefig, show

# Create simple differentiator
Ci = 1e-6
Rf = 1e3

a = Opamp()

# Connect V+ to ground.
b = a.shortcircuit(1)

# Add feedback resistor.
c = b.bridge(R(Rf))

# Add input capacitor to V-.
d = c.prepend(Series(C(Ci)))

print d.Vtransfer


f = np.logspace(1, 8, 1000)

fig = figure()
ax = fig.add_subplot(111)
Zf = d.Vtransfer.frequency_response(f)
ax.loglog(f, abs(Zf))
ax.grid(True)

show()

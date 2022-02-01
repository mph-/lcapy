import numpy as np
from matplotlib.pyplot import subplots, style, savefig, show

from lcapy import *

h = exp(-t) * u(t)
h1 = h(n, method='impulse-invariance').subs(dt, 0.1)
# Same as impulse-invariance since there are no zeros
#h2 = h(n, method='matched-Z').subs(dt, 0.1)
h3 = h(n, method='bilinear').subs(dt, 0.1)
h4 = h(n, method='forward-euler').subs(dt, 0.1)
h5 = h(n, method='backward-euler').subs(dt, 0.1)

fig, axes = subplots(1, figsize=(6, 3.5))

h1.plot((-5, 15), axes=axes, label='impulse-invariance', color='C0')
#h2.plot((-5, 15), axes=axes, label='matched-Z', color='C1')
h3.plot((-5, 15), axes=axes, label='bilinear', color='C2')
h4.plot((-5, 15), axes=axes, label='forward-Euler', color='C3')
h5.plot((-5, 15), axes=axes, label='backward-Euler', color='C4')
axes.set_ylim(0, 1.5)

axes.legend()

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')

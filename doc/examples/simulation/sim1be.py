from lcapy import Circuit
from numpy import linspace
from matplotlib.pyplot import savefig

cct = Circuit('VRL1.sch')

tv = linspace(0, 1, 100)

results = cct.sim(tv, integrator='backward-euler')

ax = cct.R1.v.plot(tv, label='analytic')
ax.plot(tv, results.R1.v, label='simulated')
ax.legend()

savefig('sim1be.png')

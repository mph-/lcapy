from lcapy import Circuit
from numpy import linspace
from matplotlib.pyplot import savefig

cct = Circuit('VRC1.sch')

tv = linspace(0, 1, 101)

results = cct.sim(tv)

ax = cct.R1.v.plot(tv, label='analytic')
ax.plot(tv, results.R1.v, label='simulated')
ax.legend()

savefig('sim2.png')

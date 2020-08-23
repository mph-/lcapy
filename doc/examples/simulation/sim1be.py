from lcapy import *
import numpy as np

a = Circuit("""
V1 1 0 step 10; down
R1 1 2 5; right
L1 2 0_2 2; down
W 0 0_2; right""")

tv = np.linspace(0, 1, 100)

sim = a.sim
results = sim(tv, integrator='backward-euler')

ax = a.R1.v.plot(tv)
ax.plot(tv, results.R1.v)

from lcapy import *
import numpy as np

a = Circuit("""
V1 1 0 step 10; down
R1 1 2 5; right
C1 2 3 0.01; right
L1 3 0_3; down
W 0 0_3; right""")

tv = np.linspace(0, 1, 100)

sim = a.sim
results = sim(tv)

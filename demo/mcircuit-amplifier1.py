from msignal.mcircuit import *
import numpy as np
from matplotlib.pyplot import figure, savefig, show

a = VoltageAmplifier()
print a.Vtransfer


a = VoltageAmplifier(10, 1e-8, 1e-10, 1e-10)
print a.Vgain(1, 2)
print a.Vgain(2, 1)



show()

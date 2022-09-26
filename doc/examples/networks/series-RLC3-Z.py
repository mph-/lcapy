from lcapy import *
from numpy import logspace
from matplotlib.pyplot import savefig

N = R(10) + C(1e-4) + L(1e-3)

vf = logspace(0, 5, 400)
N.Z(j2pif).magnitude.plot(vf)

savefig('series-RLC3-Z.png')

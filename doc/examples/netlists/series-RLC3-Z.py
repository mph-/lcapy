from lcapy import *
from numpy import logspace
from matplotlib.pyplot import savefig, show

N = R(10) + C(1e-4) + L(1e-3)

vf = logspace(0, 5, 400)
Z = N.Z.frequency_response().magnitude.plot(vf)

show()

savefig('series-RLC3-Z.png')

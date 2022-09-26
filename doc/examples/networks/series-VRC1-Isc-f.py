from lcapy import *
from numpy import logspace
from matplotlib.pyplot import savefig

N = Vstep(10) + R(10) + C(1e-4)

vf = logspace(0, 5, 400)
N.Isc(f).plot(vf, log_scale=True)

savefig('series-VRC1-Isc-f.png')

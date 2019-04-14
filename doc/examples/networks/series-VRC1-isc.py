from lcapy import *
from numpy import linspace
from matplotlib.pyplot import savefig

N = Vstep(20) + R(10) + C(1e-4, 0)

vt = linspace(0, 0.01, 1000)
N.Isc(t).plot(vt)

savefig('series-VRC1-isc.png')

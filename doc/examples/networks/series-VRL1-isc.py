from lcapy import *
from numpy import linspace
from matplotlib.pyplot import savefig

N = Vstep(20) + R(10) + L(1e-2, 0)

vt = linspace(0, 0.01, 1000)
N.Isc(t).plot(vt)

savefig('series-VRL1-isc.png')

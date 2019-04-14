from lcapy import Vstep, R, L, C, t
from matplotlib.pyplot import savefig
from numpy import linspace

a = Vstep(10) + R(0.1) + C(0.4) + L(0.2, 0)

vt = linspace(0, 10, 1000)
a.Isc(t).plot(vt)

savefig('series-VRLC1-isc.png')

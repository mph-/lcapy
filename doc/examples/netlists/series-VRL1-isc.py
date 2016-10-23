from lcapy import *
from numpy import linspace
from matplotlib.pyplot import savefig, show

N = Vstep(20) + R(10) + L(1e-2, 0)

tv = linspace(0, 0.01, 1000)
N.Isc.transient_response().plot(tv)

show()

savefig('series-VRL1-isc.png')

from lcapy import *
from numpy import linspace
from matplotlib.pyplot import savefig, show

N = Vdc(20) + R(10) + C(1e-4, 0)

tv = linspace(0, 0.01, 1000)
N.Isc.transient_response().plot(tv)

show()

savefig('series-VRC1-isc.png')

from matplotlib.pyplot import savefig
from numpy import linspace
from lcapy import *

vt = linspace(-5, 5, 200)
cos(2 * t).plot(vt)

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


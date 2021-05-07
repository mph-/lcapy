from matplotlib.pyplot import savefig
from numpy import linspace
from lcapy import *

vt = linspace(-5, 5, 200)
axes = cos(2 * t).plot(vt, label='cos')
sin(2 * t).plot(vt, label='sin', axes=axes)
axes.legend()

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


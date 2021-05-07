from matplotlib.pyplot import savefig, subplots
from numpy import linspace
from lcapy import *

figs, axes = subplots(1)
cos(2 * t).plot(axes=axes, label='cos')
sin(2 * t).plot(axes=axes, label='sin')
axes.legend()

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


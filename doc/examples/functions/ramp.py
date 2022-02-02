from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
ramp(t).plot((-2, 2), title='ramp(t)')
savefig(__file__.replace('.py', '.png'))

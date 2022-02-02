from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
rampstep(t).plot((-2, 2), title='rampstep(t)')
savefig(__file__.replace('.py', '.png'))

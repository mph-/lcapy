from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
u(t).plot((-2, 2), title='H(t)')
savefig(__file__.replace('.py', '.png'))


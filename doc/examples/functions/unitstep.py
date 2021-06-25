from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
u(n).plot((-5, 5), title='u(n)')
savefig(__file__.replace('.py', '.png'))


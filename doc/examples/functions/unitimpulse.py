from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
delta(n).plot((-5, 5), title='delta(n)')
savefig(__file__.replace('.py', '.png'))


from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
sign(n).plot((-5, 5), title='sign(n)')
savefig(__file__.replace('.py', '.png'))


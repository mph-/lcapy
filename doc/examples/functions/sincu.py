from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
sincu(t).plot((-8, 8), title='sincu(t)')
savefig(__file__.replace('.py', '.png'))


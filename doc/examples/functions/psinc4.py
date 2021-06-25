from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
psinc(4, t).plot((-2, 2), title='psinc(4, t)')
savefig(__file__.replace('.py', '.png'))


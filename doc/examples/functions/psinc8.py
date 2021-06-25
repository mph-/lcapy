from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
psinc(8, t).plot((-2, 2), title='psinc(8, t)')
savefig(__file__.replace('.py', '.png'))


from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
sincn(t).plot((-8, 8), title='sincn(t)')
savefig(__file__.replace('.py', '.png'))


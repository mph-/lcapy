from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
tri(t - 1).plot((-2, 2), title='tri(t - 1)')
savefig(__file__.replace('.py', '.png'))

from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
tri(t).plot((-2, 2), title='tri(t)')
savefig(__file__.replace('.py', '.png'))


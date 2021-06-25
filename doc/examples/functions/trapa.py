from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
trap(t, 0.5).plot((-2, 2), title='trap(t, 0.5)')
savefig(__file__.replace('.py', '.png'))


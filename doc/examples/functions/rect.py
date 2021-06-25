from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
rect(t).plot((-2, 2), title='rect(t)')
savefig(__file__.replace('.py', '.png'))


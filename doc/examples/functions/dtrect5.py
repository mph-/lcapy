from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
rect(n / 5).plot((-5, 5), title='rect(n / 5)')
savefig(__file__.replace('.py', '.png'))


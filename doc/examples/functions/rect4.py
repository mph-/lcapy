from matplotlib.pyplot import savefig, style
from lcapy import *

style.use('function.mplstyle')
rect(n / 4).plot((-5, 5), title='rect(n / 4)')
savefig(__file__.replace('.py', '.png'))


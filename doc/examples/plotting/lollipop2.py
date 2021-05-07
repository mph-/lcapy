from matplotlib.pyplot import savefig
from lcapy import *

seq((0, 0, 1, 2, 3, 0, 1, 2, 3, 0)).plot()

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


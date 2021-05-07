from matplotlib.pyplot import savefig
from lcapy import *

phasor(1 +j).plot()

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


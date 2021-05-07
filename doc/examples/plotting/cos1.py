from matplotlib.pyplot import savefig
from lcapy import *

cos(2 * t).plot()

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


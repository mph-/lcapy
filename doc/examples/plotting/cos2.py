from matplotlib.pyplot import savefig
from lcapy import *

cos(2 * t).plot((-5, 5))

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


from matplotlib.pyplot import savefig
from lcapy import *

nexpr(1).DTFT().remove_images((-5, 5)).subs(dt, 1).plot((-5, 5))

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')


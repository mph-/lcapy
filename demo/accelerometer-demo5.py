from lcapy import V, R, L, C, IdealTransformer, Shunt, Series, TSection
import numpy as np
from matplotlib.pyplot import figure, savefig, show, rcParams

rcParams.update ({'font.size': 18})
rcParams.update ({'legend.fontsize': 14})



a1 = TSection(L('m_2') | C('S_2'), L('m_1'), C('S_m') + R('R_m')).chain(IdealTransformer('alpha'))

Av1 = a1.Vgain12

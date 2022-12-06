from matplotlib.pyplot import subplots, savefig
from numpy import arange
from numpy.random import randn
from lcapy import expr

e = expr('a * exp(-t  / tau) * u(t)')
tv = arange(1, 100)
vv = e.subs({'a': 1, 'tau': 10}).evaluate(tv) + randn(len(tv)) * 0.05
results = e.estimate(tv, vv, ranges={'a': (0, 10), 'tau': (1, 20)})

vv2 = e.subs(results.params).evaluate(tv)
fig, axes = subplots(1, figsize=(6, 3))
axes.plot(tv, vv, '.')
axes.plot(tv, vv2)
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')

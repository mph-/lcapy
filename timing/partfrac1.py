import time
from lcapy import *

funcs = [1 / s, 1 / s**2, 1 / (s + 3), 1 / (s + 3)**2, (s + 3) / (s + 4),
         1 / (s + 3)**2 / (s + 4), 1 / (s + 3)**3 / (s + 4),
         1 / (s + 3) / (s + 4) / (s + 5), (s + 6) / (s + 3) / (s + 4) / (s + 5),
         1 / (s + 3)**2 / (s + 4)**2,  1 / (s + 3)**3 / (s + 4)**2,
         s / (s + 3)**2 / (s + 4), s / (s + 3)**3 / (s + 4)]

Ntrials = 10
methods = ('ec', 'sub')
times = {}

for func in funcs:
    ans1 = func.partfrac(method='ec')
    ans2 = func.partfrac(method='sub')
    if ans1 != func:
        print('Wrong answer for eq: ', func)
    if ans2 != func:
        print('Wrong answer for sub: ', func)


for method in methods:
    times[method] = []

    for func in funcs:
        start = time.perf_counter()
        for i in range(Ntrials):
            func.partfrac(method=method)
        stop = time.perf_counter()
        times[method].append((stop - start) / Ntrials)

import numpy as np
from matplotlib.pyplot import subplots, style, savefig, show

index = np.arange(len(funcs))

fig, axes = subplots(1)
axes.bar(index, times['ec'], 0.35, label='ec')
axes.bar(index+0.35, times['sub'], 0.35, label='subs')
axes.legend()
axes.set_ylabel('Time (s)')

show()

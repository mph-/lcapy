from lcapy import TP

tp1 = TP(l='Two-port 1', fill='blue')
tp2 = TP(l='Two-port 2', fill='blue')

n = tp1.parallel(tp2)

n.draw(__file__.replace('.py', '.png'))

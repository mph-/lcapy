from lcapy import TP

tp1 = TP(l='Two-port 1', fill='blue')
tp2 = TP(l='Two-port 2', fill='blue')
tp = tp1.series(tp2)

tp.draw(__file__.replace('.py', '.png'))

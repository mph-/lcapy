from lcapy import TP

tp1 = TP(label='Two-port 1')
tp2 = TP(label='Two-port 2')

n = tp1.parallel(tp2)

n.draw(__file__.replace('.py', '.png'))

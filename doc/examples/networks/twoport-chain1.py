from lcapy import TP

tp1 = TP(label='Two-port 1')
tp2 = TP(label='Two-port 2')
tp = tp1.chain(tp2)

tp.draw(__file__.replace('.py', '.png'))

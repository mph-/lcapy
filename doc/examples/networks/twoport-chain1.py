from lcapy import GenericTwoPort

tp1 = GenericTwoPort(label='Two-port 1')
tp2 = GenericTwoPort(label='Two-port 2')

n = tp1.chain(tp2)

n.draw(__file__.replace('.py', '.png'))

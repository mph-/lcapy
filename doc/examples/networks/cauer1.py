from lcapy import impedance

Z = impedance('(4*s**2 + 3 * s + 1 / 6) / (s**2 + 2 * s / 3)')

n = Z.network('cauerI')

n.draw(__file__.replace('py', 'png'), layout='ladder')

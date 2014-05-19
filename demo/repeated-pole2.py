from lcapy import s

H = (3 * s + 2) / ((s + 1)**2)

H.pprintans('H(s)')

H.partfrac().pprintans('H(s)')

H.inverse_laplace().pprintans('h(t)')

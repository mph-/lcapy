from lcapy import s

H = (2 * s**2 + 3 * s + 2) / (s + 4)

H.pprintans('H(s)')

H.partfrac().pprintans('H(s)')

H.inverse_laplace().pprintans('h(t)')


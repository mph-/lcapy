from lcapy import s

H = (s**2 + 4) / (s**2 + 3 * s + 2)

H.pprintans('H(s)')

H.partfrac().pprintans('H(s)')

H.inverse_laplace().pprintans('h(t)')


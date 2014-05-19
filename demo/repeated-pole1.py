from lcapy import s

H = 2 / (s**3 + 12 * s**2 + 36 * s)

H.pprintans('H(s)')

H.partfrac().pprintans('H(s)')

H.inverse_laplace().pprintans('h(t)')









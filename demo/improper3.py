from lcapy import s

H = (2 * s**3 + 8 * s**2 + 2*s + 4) / (s**2 + 5 * s + 6)

H.pprintans('H(s)')

H.partfrac().pprintans('H(s)')

H.inverse_laplace().pprintans('h(t)')





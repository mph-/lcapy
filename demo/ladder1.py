from mcircuit import R, L, C

Za = R(1)
Zb = R(1)
Zc = R(1)
Zd = R(1)
Ze = R(1)
Zf = R(1)

H1 = Zf.ladder(Ze)

H2 = Zf.ladder(Ze, Zd, Zc)

H3 = Zf.ladder(Ze, Zd, Zc, Zb, Za)

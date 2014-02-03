from msignal.mcircuit import R, L, C

Za = R(1)
Zb = R(1)
Zc = R(1)
Zd = R(1)
Ze = R(1)
Zf = R(1)

H1 = Zf.chain_divider(Ze)

H2 = Zf.chain_divider(Ze, Zd, Zc)

H3 = Zf.chain_divider(Ze, Zd, Zc, Zb, Za)

U1 dff; right, pinnames=all, pinlabels=all, pinnodes=all, l={dff}
U2 jkff; right, pinnames=all, pinlabels=all, pinnodes=all, l={jkff}
U3 rslatch; right, pinnames=all, pinlabels=all, pinnodes=all, l={rslatch}
U4 chip3131; right, pinlabels={l1=D,l2=>,r1=Q,r3=Q$\slash$}, pinnodes=all, l={dff}
O U1.mid U2.mid; right=2
O U2.mid U3.mid; right=2
O U3.mid U4.mid; right=3

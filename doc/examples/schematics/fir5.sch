W x 1; right
W 1 U1.in; right=0.5, endarrow=tri
U1 box .in .out; right=0.5, l=$z^{-1}$
U2 circle4 .in .mul .out .unused; down=0.5, l=$\times$
U3 circle4 .in .mul .out .unused; down=0.5, l=$\times$
U4 circle4 .add1 .unused .out .add2; right=0.5, l=$+$
W U1.out 2; right=0.5
W 3 U4.add1; right=1.25, endarrow=tri
W 1 U2.in; down
W U2.out 3; down
W a0 U2.mul; right=0.5, endarrow=tri
W 2 U3.in; down, endarrow=tri
W U3.out U4.add2; down, endarrow=tri
W a1 U3.mul; right=0.5, endarrow=tri
W U4.out y; right=0.5, endarrow=tri
# Align multipliers
O U2.mul U3.mul; right
; draw_nodes=false, label_nodes=alpha

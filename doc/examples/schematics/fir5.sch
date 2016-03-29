W 1 U1.in; right=0.5, endarrow=tri
U1 box .in .out; right=0.5, l=$z^{-1}$
W U1.out 2; right=0.5
W 3 U3.add1; right=1.25, endarrow=tri
U3 circle4 .add1 .unused .out .add2; right=0.5, l=$+$
U4 circle4 .in .mul .out .unused; down=0.5, l=$\times$
W 2 U4.in; down, endarrow=tri
W U4.out U3.add2; down, endarrow=tri
W 4 U4.mul; right=0.5, endarrow=tri
; draw_nodes=connections




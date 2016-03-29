W 1 2; right
W 2 U1.in; right=0.5, endarrow=tri
U1 box .in .out; right=0.5, l=$z^{-1}$
W U1.out 4; right=0.5
W 2 U2.in; down=0.5
U2 circle4 .in .mul .out .unused; down=0.5, l=$\times$
W 5 U2.mul; right=0.5, endarrow=tri
W U2.out 3; down
W 3 U3.add1; right=1.25, endarrow=tri
U3 circle4 .add1 .unused .out .add2; right=0.5, l=$+$



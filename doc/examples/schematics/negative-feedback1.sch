W 1 SP1._1; right=0.5, endarrow=tri
SP1 pm ._1 ._2 ._3; right=1, l={}
W SP1._3 U1._IN; endarrow=tri
#U1 box ._IN ._OUT; right=1.5, aspect=1.5, l=Open-loop gain (A)
U1 box ._IN ._OUT; right, aspect=1.5, l=Open-loop gain (A)
W U1._OUT 2; right
W 2 3; down
U2 box ._OUT ._IN; right=1.5, aspect=1.5, l=Attenuator $\beta$
W 3 U2._IN; left, endarrow=tri
W U2._OUT 4; left
W 4 SP1._2; up, endarrow=tri
W 2 5; right=0.5, endarrow=tri
; draw_nodes=false, label_nodes=false

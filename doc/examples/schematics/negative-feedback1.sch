W 1 SJ1._2; right=0.5, endarrow=tri
SJ1 _1 _2 _3 _4; right, l2=$+$, l3=$-$, l={}
W SJ1._1 U1._IN; endarrow=tri
U1 block _IN _OUT; right, l=Open-loop gain (A)
W U1._OUT 2; right
W 2 3; down
U2 block _IN _OUT; left, l=Attenuator $\beta$
W 3 U2._IN; left, endarrow=tri
W U2._OUT 4; left
W 4 SJ1._3; up, endarrow=tri
W 2 5; right=0.5, endarrow=tri
; draw_nodes=connections, label_nodes=false


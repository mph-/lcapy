W 1 SP1.1; right=0.5, endarrow=tri, l=V{in}
SP1 pm .1 .2 .3; right, l={}
W SP1.3 TR1.IN; endarrow=tri, l=Vd
TR1 .IN .OUT A; right=1.5, aspect=0.67, l=Open-loop gain (A)
W TR1.OUT 2; right, l=V{out}
W 2 3; down
TR2 .OUT .IN B; right=1.5, aspect=0.67, l=Attenuator $\beta$
W 3 TR2.IN; left, endarrow=tri
W TR2.OUT 4; left
W 4 SP1.2; up, endarrow=tri
W 2 5; right=0.5, endarrow=tri
; draw_nodes=false, label_nodes=false

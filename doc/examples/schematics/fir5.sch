W x 1; right
W 1 S1.w; right=0.5, endarrow=tri
S1 box; right=0.5, l=$z^{-1}$
S2 circle; right=0.5, l=$\times$
S3 circle; right=0.5, l=$\times$
S4 circle; right=0.5, l=$+$
W S1.e 2; right=0.5
W 3 S4.w; right=1.25, endarrow=tri
W 1 S2.n; down
W S2.s 3; down
W a0 S2.w; right=0.5, endarrow=tri
W 2 S3.n; down, endarrow=tri
W S3.s S4.n; down, endarrow=tri
W a1 S3.w; right=0.5, endarrow=tri
W S4.e y; right=0.5, endarrow=tri
# Align multipliers
O S2.w S3.w; right
; draw_nodes=false, label_nodes=alpha
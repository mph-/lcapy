TP 3 4 1 2; right, shape=cloud, l=Independent source free LTI network, fill=blue!10
W 3 m; right=0.5, i<=I_{mn}
W 4 n; right=0.5, i=I_{mn}
W i 1; right=0.5, i=I_{ij}
W j 2; right=0.5, i<=I_{ij}
A i; l=i, anchor=se
A j; l=j, anchor=ne
A m; l=m, anchor=sw
A n; l=n, anchor=nw
P1 i j; down, v=V_{ij}
P2 m n; down, v=V_{mn}
; label_style=value, draw_nodes=connected, label_nodes=none

TP 3 4 1 2; right, shape=cloud, l=Independent source free network, fill=blue!10
W 3 m; right=0.5
W 4 n; right=0.5
W i 1; right=0.5
W j 2; right=0.5
A i; l=i, anchor=se
A j; l=j, anchor=ne
A m; l=m, anchor=sw
A n; l=n, anchor=nw
I1 i j Iij; down
P1 i j; down
P2 m n; down, v=V_{mn}
; label_style=value, draw_nodes=connected, label_nodes=none

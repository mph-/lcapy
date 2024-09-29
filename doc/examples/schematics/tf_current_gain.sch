TP 3 4 1 2; right, shape=cloud, l=Independent source free LTI network, fill=blue!10
W 3 m; right=0.5
W 4 n; right=0.5
W i 1; right=0.5
W j 2; right=0.5
A i; l=i, anchor=se
A j; l=j, anchor=ne
A m; l=m, anchor=sw
A n; l=n, anchor=nw
I1 j i Iij; up
P1 i j; down
P2 m n; down
W m 5; down=0.3
W 5 6; down=0.3, i<=I_{mn}
W 6 n; down=0.3
; label_style=value, draw_nodes=connected, label_nodes=none

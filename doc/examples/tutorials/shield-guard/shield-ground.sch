Vs 14 12 ac; down
Rs 14 1; right
Cable1; right=4, dashed, kind=coax, l=
W 1 Cable1.in; right=0.5
W Cable1.out 2; right=0.5
W Cable1.ognd 10; down=0.5
Cc Cable1.mid Cable1.b; down=0.2, dashed, scale=0.6
W 2 11; right=0.75
W 7 10; up=0.5
E2 15 0 opamp 11 17 A_1; right, scale=0.5
W 17 18; down
W 12 7; right
W 7 18; right
W 18 0; down=0.2, sground
Rin 11 17; down
; label_nodes=none, draw_nodes=connections


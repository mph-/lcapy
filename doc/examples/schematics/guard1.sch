Cable1; right=4, dashed, kind=coax
W 1 Cable1.in; right=0.5
W Cable1.out 2; right=0.5
# Provide electrical connection
W Cable1.in Cable1.out; free, invisible
W Cable1.ognd 10; down=0.5
Cc Cable1.mid Cable1.b; down=0.2, dashed, scale=0.6
W 2 3; right=1.5
W 3 11; right=0.5
W 3 4; down=0.5
W 5 6; down=0.5
W 6 7; left
W 7 10; up=0.5
R 10 8; right
E 8 0 opamp 4 5 A; left=0.5, mirror, scale=0.5
; label_nodes=none, draw_nodes=connections


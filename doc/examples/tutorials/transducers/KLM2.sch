Cable1; right=2, kind=coax, l=Z_0
W 1 Cable1.in; right=0.5, i=U_b
W Cable1.out 2; right=0.5
W Cable1.ignd 1_0; down=0.5
W 14_0 1_0; right=0.5
Pb 1 14_0; down, v_=F_b

Cable2; right=2, kind=coax, l=Z_0
W 2 Cable2.in; right=0.5
W Cable2.out 3; right=0.5, i<=U_f
W Cable2.ognd 3_0; down=0.5
W 3_0 13_0; right=0.5
Pf 3 13_0; down, v^=F_f

W 2 4; down=1.5
W Cable2.ignd 5; down=0.1
W Cable1.ognd 6; down=0.1
W 6 5; right

TF1 4 7 8 9 phi; right
W 7 10; right=0.5
W 5 10; down

Pe 11 0; down, v_=V
C0 11 12; right, i>^=I
Z 12 8 {j * X_1}; right
W 0 9; right

W 9 7; right, ignore

; draw_nodes=connections, label_nodes=none, label_ids=none

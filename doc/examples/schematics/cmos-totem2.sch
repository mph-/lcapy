M1 1 2 3 pmos; right, l=Pull-up (PMOS)
M2 1 4 5 nmos; right, l=Pull-down (NMOS)
W 3 6; up=0.1, sground, l=V_{DD}
W 5 7; down=0.1, sground, l=V_{SS}
W 1 11; right=2
W 11 PIN; right, output, l=PIN
D1 11 12; up
W 12 12_1; up=0.1, implicit, l=V_{DD}
D2 13 11; up
W 13 13_1; down=0.1, implicit, l=V_{SS}

;draw_nodes=connections, label_nodes=alpha, label_ids=False

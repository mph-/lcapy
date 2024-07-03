# Created by lcapy-tk V0.95.dev0
; nodes={1@(2.5, 8), 2@(2.5, 6.5), 3@(2.5, 5), 4@(4.5, 8), 5@(4.5, 6.5), 6@(4.5, 5), 7@(6.5, 8), 8@(6.5, 6.5), 9@(6.5, 5), 10@(8.5, 8), 11@(8.5, 6.5), 12@(8.5, 5), 13@(0.5, 8), 14@(0.5, 5)}
R1 1 2; down=0.75, scale=0.75
C1 2 3; down=0.75, scale=0.75
R2 4 5; down=0.75, scale=0.75
C2 5 6; down=0.75, scale=0.75
R3 7 8; down=0.75, scale=0.75
C3 8 9; down=0.75, scale=0.75
R4 10 11; down=0.75, scale=0.75
C4 11 12; down=0.75, scale=0.75
W1 13 1; right
W2 1 4; right
W3 4 7; right
W4 7 10; right
W5 14 3; right
W6 3 6; right
W7 6 9; right
W8 9 12; right
; draw_nodes=connections, label_nodes=none, style=american, voltage_dir=RP, label_style=value

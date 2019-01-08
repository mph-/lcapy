V 1 0 VDD; down, size=1.2
R 1 2 Roh; right, size=1.5, l=R_{oh}(I_{oh}), i=I_{oh}
SW 2 3 no; right
W 3 4; right, size=0.5
C 4 5 C {-VDD / 2}; right, v=-0.5 V_{DD}
V 5 6 {0.5 * VDD}; down
W 0 7; right
W 7 6; right
P 4 7; down, v=0
;draw_nodes=connections, label_nodes=False, label_ids=False
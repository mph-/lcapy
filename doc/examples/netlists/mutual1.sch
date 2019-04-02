V 4 5 {u(t)}; down
W 4 1; right
W 5 0; right
L1 1 0; down
L2 2 3; down=1.5
# M = k * sqrt(L1 * L2)
K1 L1 L2 {M / sqrt(L1 * L2)}; size=1.5
W 2 6; right
W 3 7; right
R 6 7; down
# Need wire or high value R to link circuits
W 0 3; right
;label_ids=false

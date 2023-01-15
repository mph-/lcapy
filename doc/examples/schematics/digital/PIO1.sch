U1 buffer; right, l={}, fill=blue!50
U2 buffer; right, fliplr, l={}, fill=red!50
W U2.in 1; right=0.5
W 1 3; down=0.5
W 3 pin; bidir, l=PA0, fill=purple!50
W 2 3; up=0.5
W U1.out 2; right=0.5
W 2i U1.in; right=0.5
W 1o U2.out; right=0.5
W U1.vss 4; down=0.5
; label_nodes=none, draw_nodes=connections

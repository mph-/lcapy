U3 dff; right=1, fliplr, pinlabels={q,d,clk}, l={}, fill=red!50
W U3.d U1.q; right=0.5
U1 dff; right=1, fliplr, pinlabels={q,d,clk}, l={}, fill=red!50
U2 buffer; right, fliplr, l={}, fill=red!50
W U2.out U1.d; left=0.5
W U2.in 3; right=0.5, bidir, l=PA0, fill=purple!50
W U1.clk 4; right=0.5
W 4 5; down=0.75
W U3.clk 6; right=0.25
W 6 7; down=0.75
W 7 5; right=1.5
W 5 8; right=0.5
A 8; l={PIO clock}, anchor=south west
; label_nodes=none, draw_nodes=connections


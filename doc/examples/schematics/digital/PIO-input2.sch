U1 dff; right=1, fliplr, pinlabels={q,d,clk}, l={}, fill=red!50
U2 buffer; right, fliplr, l={}, fill=red!50
W U2.out U1.d; left=0.5
W U2.in 3; right=0.5, bidir, l=PA0, fill=purple!50
W U1.clk 4; right=0.5
W 4 5; down=0.5, input
A 5; l={PIO clock}, anchor=south west
; label_nodes=none, draw_nodes=connections

U1 dff; right=1, pinlabels={q,d,clk}, l={PORT}, fill=blue!50
U2 buffer; right, l={}, fill=blue!50
W U1.q U2.in; right=0.5
W U2.out 3; right=0.5, bidir, l=PA0, fill=purple!50

O U1.q U3.q; down=1.5
U3 dff; right=1, pinlabels={q,d,clk}, l={DD}, fill=purple!50
W U3.q U2.vss; steps=-|, free
; label_nodes=none, draw_nodes=connections



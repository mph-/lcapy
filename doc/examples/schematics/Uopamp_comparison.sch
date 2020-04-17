U1 opamp; right=1, pinnodes=all, pinnames=all, l=opamp
U2 fdopamp; right=1, pinnodes=all, pinnames=all, l=fdopamp
U3 inamp; right=1, pinnodes=all, pinnames=all, l=inamp
U4 diffamp; right=1, pinnodes=all, pinnames=all, l=diffamp, pinlabels={in+,in-}
O U1.mid U2.mid; right=2
O U2.mid U3.mid; right=2
O U3.mid U4.mid; right=2
; help_lines=0.2

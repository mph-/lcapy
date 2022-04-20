E 1 2 inamp 3 4 5 6 Ad Ac Rf; right, l=
W 5_1 5; right=0.5
W 6_1 6; right=0.5
Rg 5_1 6_1; down=0.5, scale=0.5
W 2 0; down=0.1, implicit, l={0\,\mathrm{V}}
Vs1 3_2 0_3; down
W 0_3 0; down=0.1, implicit, l={0\,\mathrm{V}}
W 3_2 3; right
Vs2 4 0_4; down
W 0_4 0; down=0.1, implicit, l={0\,\mathrm{V}}
; draw_nodes=connected

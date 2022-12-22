Q1 1 2 3 npn Q1; up, l=npn
Q2 4 5 3 pnp Q2; up, l=pnp
J1 4 6 7 njf J1; up, l=njf
J2 8 9 7 pjf J2; up, l=pjf
M1 8 10 11 nmos M1; up, l=nmos
M2 12 13 11 pmos M2; up, l=pmos
M3 12 14 15 nmos M3; up, l=nmosd
M4 16 17 15 pmos M4; up, l=pmosd

# Hack to include labels in bounding box
O 7 18; up=0.8

O 1 19; down=1.5

M5 19 20 21  M5; up, kind=nfetd, l=nfetd
M6 22 23 21  M6; up, kind=pfetd, l=pfetd
M7 22 24 25  M7; up, kind=nfet, l=nfet
M8 26 27 25  M8; up, kind=pfet, l=pfet
M9 26 28 29  M9; up, kind=nigfetd, l=nigfetd
M10 30 31 29  M10; up, kind=pigfetd, l=pigfetd
M11 30 32 33  M11; up, kind=nigfete, l=nigfete
M12 34 35 33  M12; up, kind=pigfete, l=pigfete
M13 34 36 37  M13; up, kind=nigfetebulk, l=nigfetebulk
M14 38 39 37  M14; up, kind=pigfetebulk, l=pigfetebulk

O 19 40; down=1.5

M15 40 41 42  M15; up=1.5, kind=nfet, l=nfet/bodydiode, bodydiode
M16 43 44 42  M16; up=1.5, kind=pfet, l=pfet/bodydiode, bodydiode
M17 43 45 46  M17; up=1.5, kind=nmosd, l=nmosd/bulk, bulk
M18 47 48 46  M18; up=1.5, kind=pmosd, l=pmosd/bulk, bulk

; label_nodes=none, draw_nodes=connections
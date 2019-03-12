U1 chip2222; right, pindefs={sda=b1,scl=b2}, pinlabels={sda=SDA,scl=SCL}, l={Master A}
W U1.sda 11; down
W U1.scl 12; down=0.5
U2 chip2222; right, pindefs={sda=b1,scl=b2}, pinlabels={sda=SDA,scl=SCL}, l={Master B}
W U2.sda 61; down
W U2.scl 62; down=0.5
U3 chip2222; right, pindefs={sda=b1,scl=b2}, pinlabels={sda=SDA,scl=SCL}, l={Slave}
W U3.sda 21; down
W U3.scl 22; down=0.5
W 11 61; right=2.5
W 12 62; right=2.5
W 61 21; right=2.5
W 62 22; right=2.5
W 21 31; right=2
W 22 32; right=1.75
W 31 31a; up=0.5
R1 31a 41; up=1.5
R2 32 42; up=1.5
#O 41 42; right
W 41 51; up=0.2, l={VDD}, implicit
W 42 52; up=0.2, l={VDD}, implicit
; draw_nodes=connections, label_nodes=none

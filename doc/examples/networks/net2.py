from lcapy import R, C
N = (R(1) + ((R(2) + R(3)) | (C(6) + C(7) + C(8)))) | (R(4) + R(5))
N.draw('net2.png', label_values=True)




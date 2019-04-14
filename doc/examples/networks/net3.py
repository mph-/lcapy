from lcapy import R, C, L
N = (R(1) + ((R(2) + R(3)) | (L(1) + L(2) + L(3)) | (C(6) + C(7) + C(8)))) | (R(4) + R(5))
N.draw('net3.png', label_values=True)




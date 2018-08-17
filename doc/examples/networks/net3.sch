from lcapy import R, C
N = (R(1) + ((R(2) + R(3)) | (L(1) + L(2) + L(3)) | (C(6) + C(7) + C(8)))) | (R(4) + R(5))
s = N.sch()
s.draw(label_values=True)




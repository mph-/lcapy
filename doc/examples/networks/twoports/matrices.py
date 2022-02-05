from lcapy import *
from lcapy.twoport import TPA, TPB, TPG, TPH, TPY, TPZ

V1, I1, V2, I2 = symbols('V_1 I_1 V_2 I_2')

A = TPA().params
B = TPB().params
G = TPG().params
H = TPH().params
Y = TPY().params
Z = TPZ().params


fooA = Eq(Matrix((V1, I1)), MatMul(A, Matrix((V2, -I2))), evaluate=False)

fooB = Eq(Matrix((V2, I2)), MatMul(B, Matrix((V1, -I1))), evaluate=False)

fooG = Eq(Matrix((I1, V2)), MatMul(G, Matrix((V1, I2))), evaluate=False)

fooH = Eq(Matrix((V1, I2)), MatMul(H, Matrix((I1, V2))), evaluate=False)

fooY = Eq(Matrix((I1, I2)), MatMul(Y, Matrix((V1, V2))), evaluate=False)

fooZ = Eq(Matrix((V1, V2)), MatMul(Z, Matrix((I1, I2))), evaluate=False)

print(fooA.latex())

print(fooB.latex())

print(fooG.latex())

print(fooH.latex())

print(fooY.latex())

print(fooZ.latex())

print(A.Zparams.latex())

print(A.Yparams.latex())

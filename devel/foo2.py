import numpy as np


Ax = np.matrix(((1, -1, 0, 0),
               (0,  1, -1, 0),
               (0,  0, 1, -1),
               (1,  0, 0, 0)))

bx = np.matrix(((0, -1, 0, 0),)).T

x = np.linalg.pinv(Ax) * bx


Ay = np.matrix(((1, -1, 0, 0),
               (0, 1, -1, 0),
               (0, 0, 1, -1),
               (1, 0, 0, 0)))

by = np.matrix(((-1, 0, 1, 0),)).T

y = np.linalg.pinv(Ay) * by

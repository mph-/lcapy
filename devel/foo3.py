import numpy as np

def null(A, eps=1e-15):
    """Return nullspace vector of matrix A."""

    # Add row of zeroes to bottom of matrix.
    A = np.vstack((A, np.zeros((1, A.shape[-1]))))

    U, s, Vh = np.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = np.compress(null_mask, Vh, axis=0)
    return null_space.T


Ax = np.matrix(((-1, 1, 0, 0, 0, -1, 0, 0, 0, 0),
                (0, -1, 1, 0, 0, 0, -1, 0, 0, 0),
                (0, 0, -1, 1, 0, 0, 0, -1, 0, 0),
                (-1, 0, 0, 0, 1, 0, 0, 0, -1, 0),
                (0, 0, 0, 1, -1, 0, 0, 0, 0, -1),
                (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)))

bx = np.matrix(((0, 0, 0, 0, 0, 0),)).T

Ainv = np.linalg.pinv(Ax)

x = np.linalg.pinv(Ax) * bx


null(Ax)

# Need constraints that s_n >= 1  and desire \sum s_n is minimised

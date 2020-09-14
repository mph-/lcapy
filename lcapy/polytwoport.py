"""This module provides preliminary support for polyphase two-ports.
For three-phase systems, these are actually six-ports.  They extend
the concept of cascade matrices (ABCD matrices) for polyphase systems.

Copyright 2020 Michael Hayes, UCECE

"""


from .matrix import Matrix
from .expr import expr, Expr
from .polyphase import polyphase_decompose_matrix, polyphase_compose_matrix

class Polytwoport(Matrix):

    @classmethod
    def series(cls, Za, Zb, Zc, Zg=0):

        # TODO: generalise for N phases
        Za, Zb, Zc, Zg = expr(Za), expr(Zb), expr(Zc), expr(Zg)
        
        return cls(((1, 0, 0, Za + Zg, Zg, Zg),
                    (0, 1, 0, Zg, Zb + Zg, Zg),
                    (0, 0, 1, Zg, Zg, Zc + Zg),
                    (0, 0, 0, 1, 0, 0),
                    (0, 0, 0, 0, 1, 0),
                    (0, 0, 0, 0, 0, 1)))

    
    @classmethod
    def shunt(cls, Yag, Ybg, Ycg):

        # TODO: generalise for N phases        
        Yag, Ybg, Ycg = expr(Yag), expr(Ybg), expr(Ycg)

        return cls(((1, 0, 0, 0, 0, 0),
                    (0, 1, 0, 0, 0, 0),
                    (0, 0, 1, 0, 0, 0),
                    (Yag, 0, 0, 1, 0, 0),
                    (0, Ybg, 0, 0, 1, 0),
                    (0, 0, Ycg, 0, 0, 1)))

    @property
    def N_phases(self):
        return self.shape[0] // 2

    @property
    def N(self):
        return self.N_phases

    @property
    def A(self):
        return self[0:self.N, 0:self.N]

    @property
    def B(self):
        return self[0:self.N, self.N:]

    @property
    def C(self):
        return self[self.N:, 0:self.N]

    @property
    def D(self):
        return self[self.N:, self.N:]

    @A.setter
    def A(self, A):
        self[0:self.N, 0:self.N] = A

    @B.setter
    def B(self, B):
        self[0:self.N, self.N:] = B

    @C.setter
    def C(self, C):
        self[self.N:, 0:self.N] = C

    @D.setter
    def D(self, D):
        self[self.N:, self.N:] = D

    def transform(self, W=None, Winv=None):

        if W is None:
            W = polyphase_decompose_matrix(self.N)
        
        if Winv is None:
            Winv = polyphase_decompose_matrix(self.N)

        new = self * 0
        
        new.A = Winv * self.A * W
        new.B = Winv * self.B * W
        new.C = Winv * self.C * W
        new.D = Winv * self.D * W

        return new.alpha_simplify()

    def alpha_simplify(self, alpha=None):

        if self.N != 3:
            return self
        
        return alpha_simplify3(self, alpha)

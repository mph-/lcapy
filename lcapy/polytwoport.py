"""This module provides preliminary support for polyphase two-ports.
For three-phase systems, these are actually six-ports.  They extend
the concept of cascade matrices (ABCD matrices) for polyphase systems.

Copyright 2020 Michael Hayes, UCECE

"""


from .matrix import Matrix
from .expr import expr, Expr
from .polyphase import polyphase_decompose_matrix, polyphase_compose_matrix
from .polyphase import alpha_simplify3
from sympy import Add

class Polytwoport(Matrix):

    @classmethod
    def series(cls, *args, Zn=0):
        """The args are series impedances, for example:
        
         X = Polytwoport.series('Za', 'Zb', 'Zc')

        """

        args = [expr(arg) for arg in args]
        N = len(args)

        Zn = expr(Zn).expr
        
        I = Matrix.eye(N)        
        
        obj = cls(Matrix.zeros(2 * N))
        obj.A = I
        obj.B = Matrix.diag(args) + Matrix.ones(N) * Zn
        obj.D = I
        return obj

    @classmethod
    def shunt(cls, *args):
        """The args are shunt admittances to the common node, for example:
        
         X = Polytwoport.shunt('Ya', 'Yb', 'Yc')

        """

        args = [expr(arg).expr for arg in args]        
        N = len(args)

        Ysum = Add(*args)        

        I = Matrix.eye(N)        
        C = Matrix.zeros(N)

        for i in range(N):
            for j in range(N):
                if i == j:
                    C[i, j] = args[i]
                else:
                    C[i, j] = - args[i] * args[j] / Ysum
        
        obj = cls(Matrix.zeros(2 * N))
        obj.A = I
        obj.C = C
        obj.D = I
        return obj

    @classmethod
    def star(cls, *args):
        """The args are shunt admittances to the common node, for example:
        
         X = Polytwoport.star('Yas', 'Ybs', 'Ycs', 'Ysg')

        """

        args = [expr(arg).expr for arg in args]
        N = len(args) - 1

        if N != 3:
            raise ValueError('Can only handle N=3')

        I = Matrix.eye(N)

        d1 = Add(*args[0:-1])        
        d = Add(*args)
        Yas, Ybs, Ycs, Ysg = args
        
        obj = cls(Matrix.zeros(2 * N))
        obj.A = I
        obj.B = I * 0
        obj.C = Matrix((( Yas * d1,  -Yas * Ybs, -Yas * Ycs),
                        (-Yas * Ybs,  Ybs * d1,  -Ybs * Ycs),
                        (-Yas * Ycs, -Ybs * Ycs,  Ycs * d1))) / d
        obj.D = I
        return obj
        
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

    def transform(self, W=None, Winv=None, simplify=True):

        if W is None:
            W = polyphase_decompose_matrix(self.N)
        
        if Winv is None:
            Winv = polyphase_decompose_matrix(self.N)

        new = self * 0
        
        new.A = Winv * self.A * W
        new.B = Winv * self.B * W
        new.C = Winv * self.C * W
        new.D = Winv * self.D * W

        if simplify:
            new = new.alpha_simplify()
        
        return new

    def alpha_simplify(self, alpha=None):

        if self.N != 3:
            return self
        
        return alpha_simplify3(self, alpha)

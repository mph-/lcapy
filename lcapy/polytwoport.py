from .matrix import Matrix
from .expr import expr, Expr
from .polyphase import polyphase_decompose_matrix, polyphase_compose_matrix

class Polytwoport(Matrix):

    @classmethod
    def series(cls, Za, Zb, Zc, Zn):

        # TODO: generalise for N phases
        Za, Zb, Zc, Zn = expr(Za), expr(Zb), expr(Zc), expr(Zn)
        
        return cls(((1, 0, 0, Za + Zn, Zn, Zn),
                    (0, 1, 0, Zn, Zb + Zn, Zn),
                    (0, 0, 1, Zn, Zn, Zc + Zn),
                    (0, 0, 0, 1, 0, 0),
                    (0, 0, 0, 0, 1, 0),
                    (0, 0, 0, 0, 0, 1)))

    
    @classmethod
    def shunt(cls, Yan, Ybn, Ycn):

        # TODO: generalise for N phases        
        Yan, Ybn, Ycn = expr(Yan), expr(Ybn), expr(Ycn)

        return cls(((1, 0, 0, 0, 0, 0),
                    (0, 1, 0, 0, 0, 0),
                    (0, 0, 1, 0, 0, 0),
                    (Yan, 0, 0, 1, 0, 0),
                    (0, Ybn, 0, 0, 1, 0),
                    (0, 0, Ycn, 0, 0, 1)))

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

    
def alpha_simplify3(self, alpha=None):    
        
    if alpha is None:
        alpha = expr('alpha')
        
    new1 = self.expand()
    new2 = new1.replace(alpha**4, alpha)
    new3 = new2.replace(alpha**3, 1)
    new4 = new3.replace(alpha**2, -1 - alpha)
    new = new4.simplify()
    
    return new
    

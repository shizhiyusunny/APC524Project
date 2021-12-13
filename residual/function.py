from fenics import *
class FunctionW:
    def __init__(self, f):
        u, c_trans, c_rot = split(f)
        self.displacement = u
        self.lm_trans = c_trans
        self.lm_rot = c_rot
        self.dim = 2
        self.unknown = f

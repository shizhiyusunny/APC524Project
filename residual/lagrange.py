from fenics import *
from . import energy

class LagrangeMultiplier(energy.EnergyFunctional):
    def __init__(self, inputs):
        self.dof = inputs['dof']
    def return_energy(self, f):
        c_trans = f.lm_trans
        c_rot = f.lm_rot
        u = f.displacement
        if self.dof == 'translation':
            constraint = dot(c_trans,u)*dx
        elif self.dof == 'rotation':
            r=Expression(('x[0]','x[1]'),degree=1)
            constraint = c_rot*(r[0]*u[1]-r[1]*u[0])*dx
        return constraint

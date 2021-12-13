from fenics import *
from . import energy

class Traction(energy.EnergyFunctional):
    def __init__(self, inputs, mesh):
        self.traction = Traction.create_traction_object(inputs['sigma'], mesh)

    def create_traction_object(trac, mesh):
        sigma = Constant(((trac['xx'], trac['xy']), (trac['xy'], trac['yy'])))
        n = FacetNormal(mesh)
        traction = sigma * n
        return traction

    def return_energy(self, f):
        u = f.displacement
        virtual_work = -dot(self.traction,u)*ds
        return virtual_work

    def update(self, val):
        pass

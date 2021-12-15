from fenics import *
from .. import energy

class Elasticity(energy.EnergyFunctional):
    def handler(inputs, material_constants):
        from . import linear
        return linear.LinearElasticity(material_constants)
    def return_energy(self, f):
        pass

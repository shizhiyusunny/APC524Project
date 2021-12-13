from fenics import *
from . import energy

class Elasticity(energy.EnergyFunctional):
    def handler(inputs, material_constants):
        return LinearElasticity(material_constants)
    def return_energy(self, f):
        pass

class LinearElasticity(Elasticity):
    def __init__(self, material_constants):
        self.mu = material_constants['mu']
        self.Lambda = material_constants['lambda']
    def return_energy(self, f):
       u = f.displacement
       d = f.dim
       epsilon = sym(grad(u))
       sigma = 2*self.mu*epsilon + self.Lambda*tr(epsilon)*Identity(d)
       elastic_energy = 1/2*inner(sigma,epsilon)*dx
       return elastic_energy

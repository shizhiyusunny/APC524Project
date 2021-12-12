from fenics import *

class EnergyFunctional:
    def handler(mesh, material_constants):

        #clean this up.
        sigma_xx = 0.2
        sigma_xy = 0
        sigma_yy = 0
        sigma_0 = Constant(((sigma_xx,sigma_xy),(sigma_xy,sigma_yy)))
        n = FacetNormal(mesh)
        #

        energy_functionals = []
        tractions = []
        energy_functionals.append(LinearElasticity(material_constants))
        traction = Traction(sigma_0 * n)
        energy_functionals.append(traction)
        tractions.append(traction)
        energy_functionals.append(LagrangeMultiplier())
        return energy_functionals, tractions
    def return_energy(self, f):
        pass

class LinearElasticity(EnergyFunctional):
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

class LagrangeMultiplier(EnergyFunctional):
    def return_energy(self, f):
        c_trans = f.lm_trans
        c_rot = f.lm_rot
        u = f.displacement
        r=Expression(('x[0]','x[1]'),degree=1)
        constraints = dot(c_trans,u)*dx + c_rot*(r[0]*u[1]-r[1]*u[0])*dx
        return constraints

class Traction(EnergyFunctional):
    def __init__(self, t):
        self.traction = t
    def return_energy(self, f):
        u = f.displacement
        virtual_work = -dot(self.traction,u)*ds
        return virtual_work

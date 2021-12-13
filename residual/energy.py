from fenics import *

class EnergyFunctional:
    def handler(mesh, inputs, material_constants):
        energy_functionals = []
        tractions = []
        for key, val in inputs.items():
            if key == 'Elasticity':
                energy_functionals.append(Elasticity.handler(val, material_constants))
                print('Created elasticity')
            elif key == 'Tractions':
                for trac in val:
                    traction = Traction(trac, mesh)
                    energy_functionals.append(traction)
                    tractions.append(traction)
                print('Created tractions')
            elif key == 'Constraints':
                for v in val:
                    if 'Lagrange Multiplier' in v:
                        energy_functionals.append(LagrangeMultiplier(v['Lagrange Multiplier']))
                print('Created constraints')
        return energy_functionals, tractions
    def return_energy(self, f):
        pass

class Elasticity(EnergyFunctional):
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

class LagrangeMultiplier(EnergyFunctional):
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

class Traction(EnergyFunctional):
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

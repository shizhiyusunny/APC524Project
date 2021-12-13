from fenics import *

class EnergyFunctional:
    def handler(mesh, inputs, material_constants):
        from . import elasticity, traction, lagrange
        energy_functionals = []
        tractions = []
        for key, val in inputs.items():
            if key == 'Elasticity':
                energy_functionals.append(elasticity.Elasticity.handler(val, material_constants))
                print('Created elasticity')
            elif key == 'Tractions':
                for trac in val:
                    traction = traction.Traction(trac, mesh)
                    energy_functionals.append(traction)
                    tractions.append(traction)
                print('Created tractions')
            elif key == 'Constraints':
                for v in val:
                    if 'Lagrange Multiplier' in v:
                        energy_functionals.append(lagrange.LagrangeMultiplier(v['Lagrange Multiplier']))
                print('Created constraints')
        return energy_functionals, tractions
    def return_energy(self, f):
        pass

from . import energy
from fenics import *
class Residual:
    def __init__(self, mesh, material_constants):
        self.energy_functionals, self.tractions = energy.EnergyFunctional.handler(mesh, material_constants)

    def calculate_residual(self, f):
        first = True
        for energy_functional in self.energy_functionals:
            if first:
                self.free_energy = energy_functional.return_energy(f)
                first = False
            else:
                self.free_energy += energy_functional.return_energy(f)
        self.unknown = f.unknown
        self.Res = derivative(self.free_energy, self.unknown)

    def solve(self):
        solve(self.Res == 0, self.unknown)

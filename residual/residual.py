from . import energy
from . import funcw
from fenics import *
class Residual:
    def builder(inputs, mesh, material_constants):
        residual = Residual(inputs, mesh, material_constants)
        func = funcw.FunctionW(inputs['displacement'], mesh, residual.lm)
        return func, residual
    def __init__(self, inputs, mesh, material_constants):
        self.energy_functionals, self.tractions, self.lm = energy.EnergyFunctional.handler(mesh, inputs, material_constants)

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

    def update_tractions(self, val):
        for t in self.tractions:
            t.update(val)

from __future__ import print_function
from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from residual import residual
from annealer import annealer
from material import material
import yaml

# Create rectangular mesh with two circular inclusions
N = 100
L = 1
R_2 = 0.05
R_3 = 0.08
domain = Rectangle(Point(-L/2,-L/2),Point(L/2,L/2))
# mark subdomains with markers 1, 2, 3
domain.set_subdomain(1, Rectangle(Point(-L/2,-L/2),Point(L/2,L/2)))
domain.set_subdomain(2, Circle(Point(0.,0.43), R_2))
domain.set_subdomain(3, Circle(Point(-0.15,0.35), R_3))
mesh = generate_mesh(domain, N)
d = mesh.topology().dim() # dimensionality of the problem
markers = MeshFunction("size_t", mesh, d , mesh.domains())

#read yaml file
with open('test.yml') as file:
    inputs = yaml.load(file, Loader=yaml.FullLoader)

material_constants = material.MaterialConstant.builder(markers, inputs['Material Constant'])
func, residual = residual.Residual.builder(inputs['Equilibrium'], mesh, material_constants)
annealer = annealer.Annealer.build(inputs['Annealer'])
annealer.stepper(residual, func)
print("Free energy:", assemble(residual.free_energy))

# export displacements
VFS = VectorFunctionSpace(mesh, 'Lagrange', 1)
disp=project(func.displacement, VFS)
disp.rename("displacements","")
fileD = File("data/tractions_displacement.pvd");
fileD << disp;

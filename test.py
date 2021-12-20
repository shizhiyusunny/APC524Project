from __future__ import print_function
from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from residual import residual
from annealer import annealer
from material import material
from mesh import meshBuilder
import yaml

#read yaml file
with open('test.yml') as file:
    inputs = yaml.load(file, Loader=yaml.FullLoader)

mesh, markers = meshBuilder.MeshBuilder.build(inputs['Mesh'])
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
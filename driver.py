from __future__ import print_function
from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from residual import residual
from annealer import annealer
from material import material
from postprocess import postProcessing
from mesh import meshgenerator
from optimization import optimizer
import yaml

#read yaml file
with open('test.yml') as file:
    inputs = yaml.load(file, Loader=yaml.FullLoader)
not_converged = True
optimizer = optimizer.Optimizer.builder(inputs)

while not_converged:
    print('here')
    meshgen = meshgenerator.MeshGenerator(inputs['Mesh'])
    mesh, markers = meshgen.create_mesh_object()
    material_constants = material.MaterialConstant.builder(markers, inputs['Material Constant'])
    func, res = residual.Residual.builder(inputs['Equilibrium'], mesh, material_constants)
    anneal = annealer.Annealer.build(inputs['Annealer'])
    processes = postProcessing.Postprocessing.build(inputs['Post-Processing'])
    anneal.stepper(res, func, processes)
    not_converged = optimizer.check(func)
    print('Free Energy:', assemble(res.free_energy))

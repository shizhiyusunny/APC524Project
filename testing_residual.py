from __future__ import print_function
from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from residual import residual
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


# define boundary subdomains
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], -L/2)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], +L/2)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], L/2)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], -L/2)

left = Left()
right = Right()
top = Top()
bottom = Bottom()

# mark boundary subdomains with markers 1, 2, 3, 4
boundaries = MeshFunction("size_t", mesh, d-1, 0)
boundaries.set_all(0)
left.mark(boundaries, 1)
right.mark(boundaries, 2)
top.mark(boundaries, 3)
bottom.mark(boundaries, 4)

# elastic constants of the matrix and two circular inclusions
E_1 = 1
E_2 = 10
E_3 = 0.1
nu_1 = 0.3
nu_2 = 0.2
nu_3 = 0.1

# define class for calculating the Young's modulus over the whole domain
class E_class(UserExpression):
    def __init__(self, **kwargs):
        self.markers = markers
        super().__init__(**kwargs)
    def eval_cell(self, value, x, ufc_cell):
        if markers[ufc_cell.index] == 1:
            value[0] = E_1
        elif markers[ufc_cell.index] == 2:
            value[0] = E_2
        else:
            value[0] = E_3
    def value_shape(self):
        return ()

# define class for calculating the Poisson's ratio over the whole domain
class nu_class(UserExpression):
    def __init__(self, **kwargs):
        self.markers = markers
        super().__init__(**kwargs)
    def eval_cell(self, value, x, ufc_cell):
        if markers[ufc_cell.index] == 1:
            value[0] = nu_1
        elif markers[ufc_cell.index] == 2:
            value[0] = nu_2
        else:
            value[0] = nu_3
    def value_shape(self):
        return ()

# functions of elastic constants on the whole domain
E = E_class(degree=1)
nu = nu_class(degree=1)
mu = E/2/(1+nu)
Lambda = E*nu/(1-nu*nu)

material_constants = {'mu':mu, 'lambda':Lambda}

#read yaml file
with open('testr.yml') as file:
    inputs = yaml.load(file, Loader=yaml.FullLoader)

func, residual = residual.Residual.builder(inputs, mesh, material_constants)
residual.calculate_residual(func)
residual.solve()

#calculate total free energy
print("Tot Free Energy = ",assemble(residual.free_energy))

# export displacements
VFS = VectorFunctionSpace(mesh, 'Lagrange', 1)
disp=project(func.displacement, VFS)
disp.rename("displacements","")
fileD = File("data/tractions_displacement.pvd");
fileD << disp;
"""
# calculate and export von Mises stress
FS = FunctionSpace(mesh, 'Lagrange', 1)
devStress = sigma(u) - (1./d)*tr(sigma(u))*Identity(d)  # deviatoric stress
von_Mises = project(sqrt(3./2*inner(devStress, devStress)), FS)
von_Mises.rename("von Mises","")
fileS = File("data/tractions_vonMises_stress.pvd");
fileS << von_Mises;

# calculate and export stress component sigma_xx
sigma_xx = project(sigma(u)[0,0], FS)
sigma_xx.rename("sigma_xx","")
fileS = File("data/tractions_sigma_xx.pvd");
fileS << sigma_xx;

# calculate and export stress component sigma_yy
sigma_yy = project(sigma(u)[1,1], FS)
sigma_yy.rename("sigma_yy","")
fileS = File("data/tractions_sigma_yy.pvd");
fileS << sigma_yy;

# calculate and export stress component sigma_xy
sigma_xy = project(sigma(u)[0,1], FS)
sigma_xy.rename("sigma_xy","")
fileS = File("data/tractions_sigma_xy.pvd");
fileS << sigma_xy;

# export Young's modulus
young = project(E, FS)
young.rename("Young's modulus","")
fileS = File("data/tractions_young.pvd");
fileS << young;
"""

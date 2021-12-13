#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from fenics import *
import numpy as np
import constant as cst #Access to constant
class Simulator:
    def __init__(self,mesh):
        self.mesh=mesh
        self.strain=None
        self.stress=None
        self.disp=None
    def setFunctionSpace(self):
        degreeElements = 1
        P1 = FiniteElement('Lagrange', self.mesh.ufl_cell(), degreeElements)
        R = FiniteElement('Real', self.mesh.ufl_cell(), 0)
        MFS = FunctionSpace(self.mesh, MixedElement([(P1*P1),(R*R),R]))
        return MFS
    def straincalc(u):
        strain=sym(grad(u))
        return strain
    def stresscalc(u):
        stress=2*cst.mu*epsilon(u) + cst.lambda*tr(epsilon(u))*Identity(cst.d) # Assume we know mu, lambda, and d from global constants
        return stress
    def dispcalc(u):
        VFS = VectorFunctionSpace(mesh, 'Lagrange', 1)
        disp=project(u, VFS)
        return disp
    def stepForward(self, Res):
        MFS=setFunctionSpace()
        f = Function(MFS)
        tf = TestFunction(MFS)
        solve(Res == 0, f, tf)
        u, c_trans, c_rot = split(f)
        self.strain=straincalc(u)
        self.stress=stresscalc(u)
        self.disp=dispcalc(u)
        return u, c_trans, c_rot, self.strain, self.stress, self.disp
        


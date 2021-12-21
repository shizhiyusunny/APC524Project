from . import postProcessing
from fenics import *

class StorePVD(postProcessing.Postprocessing):
    def __init__(self, iterative):
        self.file = File("data/displacement.pvd")
        self.iterative = iterative
    def store(self, f, t=1.):
        print("Time = ", t)
        VFS = VectorFunctionSpace(f.mesh, 'Lagrange', 1, f.dim)
        u_store=project(f.displacement, VFS)
        u_store.rename("displacement","")
        if self.iterative:
            self.file << (u_store, t)
        else:
            self.file << u_store


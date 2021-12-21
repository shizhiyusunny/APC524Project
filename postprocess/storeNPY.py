from . import postProcessing
from fenics import *

class storeNPY(postProcessing.Postprocessing):
    def store(self, f, iterative, t=1):
        u=f.displacement
        arr_u=u.compute_vertex_values()
        fileD = File("data/displacement.txt")
        fileD << arr_u

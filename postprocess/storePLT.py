from . import postProcessing
from fenics import *
import matplotlib.pyplot as plt

class StorePLT(postProcessing.Postprocessing):
    def store(self, f, t=1.):
        a=plt.figure(0)
        plot(f.mesh,linewidth=0.3)
        plot(f.displacement, mode='color')
        plt.savefig('data/displacement.png')


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from fenics import *
import numpy as np
class Visualizer:
    def VisualizeMesh(mesh):
        plot(mesh,linewidth=0.3)
        plt.show()
        return 1
    def VisualizeDisplacement(u,mesh):
        VisulaizeMesh(mesh)
        plot(u,mode='color')


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from fenics import *
import numpy as np
class Exporter:
    def storeNumpy(u):
        np.save("myFile",u)
        return u
    def storePVD(u):
        u.rename("displacements","")
        fileD = File("data/clamped_displacement.pvd");
        fileD << u;


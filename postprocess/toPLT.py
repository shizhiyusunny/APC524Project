#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from . import Postprocessing
from fenics import *
import matplotlib.pyplot as plt

class toDisPLT(Postprocessing.Postprocessing):
    def store(t,f,mesh):
        a=plt.figure(0)
        plt.plot(mesh,linewidth=0.3)
        plt.plot(f.displacement, mode='color')
        plt.savefig('figure/displacement.png')

